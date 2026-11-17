function mpc = Converte_MATPOWER(m, mapaAreas, combGCI)
%CONVERTE_MATPOWER  Converte um caso do MATPOWER para o formato do codigo.
%
%   mpc = Converte_MATPOWER(m, mapaAreas, combGCI)
%
%   m         : caso MATPOWER (com ou sem o campo genfuel)
%   mapaAreas : vetor em que mapaAreas(a) e a area de controle da area a do
%               caso original (permite agrupar areas). Se tiver um valor por
%               barra (numel = numero de barras), e a area de cada barra.
%               [] => mantem as areas do caso
%   combGCI   : combustiveis que podem ser GCI (ex.: {'ng','coal','hydro'})
%
%   Formato de saida (por barra, numeracao sequencial 1..N):
%     bus : [i tipo Pd Qd Gs Bs area Vm Va baseKV ctrl Vmax Vmin PG QG Qmin Qmax]
%     gen : uma linha por barra com geracao (unidades ligadas agregadas)
%           [barra Pg Qg Qmax Qmin Vg mBase status Pmax Pmin 0 ... 0]
%     branch: formato MATPOWER, com tap 0 substituido por 1
%   Campos adicionais:
%     mpc.numOrig : numero original da barra
%     mpc.rateA   : capacidade dos ramos [MVA]
%
%   A coluna ctrl (11) indica as barras candidatas a GCI: todas as unidades
%   ligadas da barra usam combustiveis de combGCI.



    %% ____________________________________________________________________
    %%
    %%                        1. BARRAS
    %% ____________________________________________________________________

    %% a) Renumeracao
    nb      = size(m.bus, 1);
    numOrig = m.bus(:,1);
    mapa    = zeros(max(numOrig), 1);
    mapa(numOrig) = 1:nb;


    %% b) Areas de Controle
    areaOrig = m.bus(:,7);

    if isempty(mapaAreas)
        area = areaOrig;
    elseif numel(mapaAreas) == nb
        area = mapaAreas(:);
    else
        area = mapaAreas(areaOrig);
        area = area(:);
    end


    %% c) Unidades Ligadas
    on    = m.gen(:,8) == 1;
    gbus  = mapa(m.gen(on,1));
    g     = m.gen(on,:);

    if isfield(m, 'genfuel')
        fuel = m.genfuel(on);
    else
        fuel = repmat({'ng'}, sum(on), 1);
    end

    desp = ismember(fuel, combGCI);


    %% d) Agregacao por Barra
    PG     = accumarray(gbus, g(:,2), [nb 1]);
    Qmax   = accumarray(gbus, g(:,4), [nb 1]);
    Qmin   = accumarray(gbus, g(:,5), [nb 1]);
    Pmax   = accumarray(gbus, g(:,9), [nb 1]);
    Pmin   = accumarray(gbus, g(:,10), [nb 1]);
    nUn    = accumarray(gbus, 1, [nb 1]);
    nDesp  = accumarray(gbus, double(desp), [nb 1]);
    temGer = nUn > 0;

    Vg = m.bus(:,8);

    for i = 1:numel(gbus)
        Vg(gbus(i)) = g(i,6);
    end


    %% e) Tipos de Barra
    tipo = m.bus(:,2);
    tipo(tipo == 2 & ~temGer) = 1;
    tipo(tipo == 1 & temGer)  = 2;
    tipo(tipo == 4) = 1;


    %% f) Candidatas a GCI: todas as unidades da barra sao despachaveis
    ctrl = double(temGer & nDesp == nUn & tipo == 2);


    %% g) Matriz bus no Formato do Codigo
    mpc.version = '2';
    mpc.baseMVA = m.baseMVA;
    mpc.bus = [(1:nb).', tipo, m.bus(:,3), m.bus(:,4), m.bus(:,5), m.bus(:,6), area, ...
               Vg, m.bus(:,9), m.bus(:,10), ctrl, m.bus(:,12), m.bus(:,13), ...
               PG, zeros(nb,1), Qmin, Qmax];



    %% ____________________________________________________________________
    %%
    %%                   2. GERADORES E RAMOS
    %% ____________________________________________________________________

    %% a) Geradores Agregados por Barra
    bg = find(temGer);
    ng = numel(bg);
    mpc.gen = zeros(ng, 21);
    mpc.gen(:,1)  = bg;
    mpc.gen(:,2)  = PG(bg);
    mpc.gen(:,4)  = Qmax(bg);
    mpc.gen(:,5)  = Qmin(bg);
    mpc.gen(:,6)  = Vg(bg);
    mpc.gen(:,7)  = m.baseMVA;
    mpc.gen(:,8)  = 1;
    mpc.gen(:,9)  = Pmax(bg);
    mpc.gen(:,10) = Pmin(bg);


    %% b) Ramos Ligados
    br = m.branch(m.branch(:,11) == 1, :);
    br(:,1) = mapa(br(:,1));
    br(:,2) = mapa(br(:,2));
    br(br(:,9) == 0, 9) = 1;
    mpc.branch = br(:,1:13);


    %% c) Campos Adicionais
    mpc.numOrig = numOrig;
    mpc.rateA   = mpc.branch(:,6);

    if isfield(m, 'bus_name')
        mpc.bus_name = m.bus_name;
    end
end
