function D = CIA_Dados(mpc, cfg, cen)
%CIA_DADOS  Monta a estrutura de dados usada pelo controle de intercambio.
%
%   D = CIA_Dados(mpc, cfg)
%   D = CIA_Dados(mpc, cfg, cen)
%
%   mpc: caso no formato de Converte_MATPOWER
%   cfg: configuracao (Config_CIA)
%   cen: cenario [fCarga, fEol_1, ..., fEol_m], em que fCarga multiplica as
%        cargas ativa e reativa e fEol_i soma fEol_i*cfg.WF_MW(i) a
%        geracao da barra cfg.WF(i). Padrao: [1, 0, ..., 0].
%
%   Todas as grandezas de D estao em p.u. na base do caso.



    %% ____________________________________________________________________
    %%
    %%                        1. DADOS DO CENARIO
    %% ____________________________________________________________________

    %% a) Valores Padrao
    if nargin < 3 || isempty(cen)
        cen = [1, zeros(1, numel(cfg.WF))];
    end

    baseMVA = mpc.baseMVA;
    bus     = mpc.bus;
    gen     = mpc.gen;
    branch  = mpc.branch;
    fCarga  = cen(1);
    fEol    = cen(2:end);


    %% b) Carga
    bus(:,3) = mpc.bus(:,3)*fCarga;
    bus(:,4) = mpc.bus(:,4)*fCarga;


    %% c) Variacao da Geracao Eolica
    for i = 1:numel(cfg.WF)
        bus(cfg.WF(i), 14) = mpc.bus(cfg.WF(i), 14) + fEol(i)*cfg.WF_MW(i);
    end



    %% ____________________________________________________________________
    %%
    %%                     2. ESTRUTURA PARA O CONTROLE
    %% ____________________________________________________________________

    %% a) Matriz de Admitancia Nodal (esparsa)
    FE = 1;
    FC = 1;
    YBus = matriz_YBus(baseMVA, bus, gen, branch, FE, FC);

    [De, Para, ~, ~, ~, Y, ~, SHL, ~, ~, ~, ~, ~, ~, ~, ~, ~, Pd, ~, Qd, ~, ~, N, ~] = ...
        ler_dados(baseMVA, bus, gen, branch, FE, FC);


    %% b) Topologia e Tipos de Barra
    D.N    = N;
    D.Ybus = YBus;
    D.Ysh  = full(sum(YBus, 2));
    D.Ysr  = YBus - spdiags(D.Ysh, 0, N, N);

    D.tipo = bus(:,2);
    D.ref  = find(D.tipo == 3);
    D.tipo(D.ref) = 0;

    D.Vesp = bus(:,8);
    D.area = bus(:,7);
    D.Pd   = Pd;
    D.Qd   = Qd;


    %% c) Injecoes Liquidas Especificadas
    D.P0 = bus(:,14)/baseMVA - Pd;
    D.Q0 = bus(:,15)/baseMVA - Qd;


    %% d) Areas Controladas e Interligacoes
    % Todas as areas, exceto a da barra de referencia, tem restricao de
    % intercambio. Em cada interligacao da area a, r e o terminal dentro de
    % a e t o terminal fora. NI_a e a potencia ativa recebida pela area.
    % Modelo pi do ramo com tap a no lado De: serie y/a; shunts
    % y/a^2 - y/a (De) e y - y/a (Para), mais o carregamento j*b/2.
    D.areaRef   = D.area(D.ref);
    D.areasCtrl = setdiff(unique(D.area), D.areaRef).';
    D.nA        = numel(D.areasCtrl);
    lt          = find(D.area(De) ~= D.area(Para));

    for ia = 1:D.nA
        a  = D.areasCtrl(ia);
        la = lt(D.area(De(lt)) == a | D.area(Para(lt)) == a);

        D.lig(ia).r   = zeros(numel(la), 1);
        D.lig(ia).t   = zeros(numel(la), 1);
        D.lig(ia).y   = zeros(numel(la), 1);
        D.lig(ia).ysh = zeros(numel(la), 1);

        for i = 1:numel(la)
            L  = la(i);
            aL = branch(L,9);

            if D.area(Para(L)) == a
                D.lig(ia).r(i)   = Para(L);
                D.lig(ia).t(i)   = De(L);
                D.lig(ia).ysh(i) = Y(L) - Y(L)/aL + SHL(L);
            else
                D.lig(ia).r(i)   = De(L);
                D.lig(ia).t(i)   = Para(L);
                D.lig(ia).ysh(i) = Y(L)/aL^2 - Y(L)/aL + SHL(L);
            end

            D.lig(ia).y(i) = Y(L)/aL;
        end
    end


    %% e) Todos os Ramos (limite termico)
    % Capacidade em MVA (Inf quando nao informada no caso).
    aTap           = branch(:,9);
    D.ramo.de      = De;
    D.ramo.para    = Para;
    D.ramo.y       = Y./aTap;
    D.ramo.yshDe   = Y./aTap.^2 - Y./aTap + SHL;
    D.ramo.yshPara = Y - Y./aTap + SHL;
    D.ramo.rate    = mpc.rateA(:);
    D.ramo.rate(D.ramo.rate <= 0) = Inf;


    %% f) GCIs: Geracao Programada e Limites de P
    D.gci    = union(find(bus(:,11) == 1), cfg.GCIextra(:));
    D.Pg_nom = bus(:,14)/baseMVA;
    D.Pmax   = zeros(N,1);
    D.Pmin   = zeros(N,1);

    for k = D.gci.'
        g = find(gen(:,1) == k, 1);
        D.Pmax(k) = gen(g,9)/baseMVA;
        D.Pmin(k) = gen(g,10)/baseMVA;
    end


    %% g) Limites de Potencia Reativa dos Geradores (soma por barra)
    D.Qmax    = inf(N,1);
    D.Qmin    = -inf(N,1);
    barrasGer = unique(gen(:,1));

    for k = barrasGer.'
        g = (gen(:,1) == k);
        D.Qmax(k) = sum(gen(g,4))/baseMVA;
        D.Qmin(k) = sum(gen(g,5))/baseMVA;
    end


    %% h) Intercambios Especificados (um por area controlada)
    if numel(cfg.NIesp_MW) ~= D.nA
        error('CIA_Dados: cfg.NIesp_MW deve ter %d valores (areas %s).', D.nA, mat2str(D.areasCtrl));
    end

    D.NIesp = cfg.NIesp_MW(:)/baseMVA;


    %% i) Geracao Programada da Barra de Referencia (injecao liquida)
    if isempty(cfg.PrefEsp_MW)
        D.PrefEsp = D.P0(D.ref);
    else
        D.PrefEsp = cfg.PrefEsp_MW/baseMVA - Pd(D.ref);
    end


    %% j) Base e Configuracao
    D.baseMVA = baseMVA;
    D.cfg     = cfg;
end
