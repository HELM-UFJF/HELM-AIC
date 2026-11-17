function res = CIA_Controle(D, metodo)
%CIA_CONTROLE  Laco externo do controle de intercambio multiarea com limites.
%
%   metodo = 'HELM' ou 'NR'. Ambos usam exatamente o mesmo laco externo, de
%   modo que a comparacao isola o metodo de solucao da estrategia de controle.
%
%   Cada estagio resolve o FP com as restricoes NI_a = NIesp_a das areas
%   controladas que ainda tem GCI ativo e, depois da convergencia, verifica
%   os limites na ordem D.cfg.ordemLim ('QP' ou 'PQ'):
%     Q) Limites de Q (se D.cfg.limQ): ate D.cfg.nQporEstagio barras PV/H
%        violadas (as de maior violacao) passam a ter Q fixo no limite
%        (PV -> PQ; H -> tipo 4, que mantem o ajuste de P).
%     P) Limites de P dos GCIs, avaliados com lambda(1):
%          'um'   : fixa apenas o GCI de maior violacao;
%          'todos': fixa todos os violados de uma vez (como em Luo e Liu, 2023).
%        Os fatores beta da area sao renormalizados.
%     T) Limite termico (se D.cfg.limTermico): ver CIA_VerificaTermico.
%   Qualquer mudanca reinicia o estagio pela solucao germinal (HELM) ou pela
%   partida plana ou quente (NR, conforme D.cfg.partidaNR).
%
%   Perda de controlabilidade (por area): quando uma area controlada fica sem
%   GCI ativo, sua restricao e seu lambda saem do problema e o intercambio
%   alcancado e reportado; as demais areas continuam controladas.
%   Se a area da folga fica sem GCI ativo, a folga absorve a compensacao.



    %% ____________________________________________________________________
    %%
    %%                       1. INICIALIZACAO
    %% ____________________________________________________________________

    %% a) Estado do Laco
    cfg    = D.cfg;
    N      = D.N;

    % GCI fixado no limite de P
    fixP   = false(N,1);

    % Valor de Pg fixado [p.u.]
    PfixV  = zeros(N,1);

    % Barra com Q fixado no limite
    fixQ   = false(N,1);

    % Valor de Qg fixado [p.u.]
    QfixV  = zeros(N,1);

    % Retornos PQ -> PV de cada barra
    nRetQ  = zeros(N,1);

    % [barra, estagio] dos retornos
    regRet = [];

    % [barra, limite, violacao, estagio]
    regP   = [];

    % [barra, limite, violacao, estagio]
    regQ   = [];
    Cb0    = [];

    % Ramo monitorado em cada area (0 = NI)
    termL  = zeros(1, D.nA);

    % Terminal monitorado (barra)
    termE  = zeros(1, D.nA);

    % Fluxo ativo alvo no terminal [p.u.]
    termP  = zeros(1, D.nA);

    % Ajustes/trocas por area
    nTerm  = zeros(1, D.nA);

    % [area, ramo, Plim(MW), estagio]
    regT   = [];

    % Avisos do limite termico
    infoT  = {};
    V0     = [];
    est    = 0;

    % Soma de ordens (HELM) ou iteracoes (NR)
    custo  = 0;



    %% ____________________________________________________________________
    %%
    %%                     2. LACO DE ESTAGIOS
    %% ____________________________________________________________________

    while true
        est = est + 1;

        if est > cfg.maxEst
            error('CIA_Controle: numero maximo de estagios excedido.');
        end


        %% a) Tipos de Barra e Injecoes do Estagio
        tipo = D.tipo;
        P    = D.P0;
        Q    = D.Q0;

        idxP = find(fixP);
        tipo(idxP) = 2;
        P(idxP)    = PfixV(idxP) - D.Pd(idxP);

        ativos = setdiff(D.gci(:), idxP);


        %% b) Areas Controladas Ainda Ativas
        temGCI = false(1, D.nA);

        for ia = 1:D.nA
            temGCI(ia) = any(D.area(ativos) == D.areasCtrl(ia));
        end

        areasAtivas = find(temGCI);

        % Perda de controlabilidade por area
        loc         = ~temGCI;


        %% c) Coeficientes de Participacao e Estrutura do Controle
        [Cb, temRef] = CIA_Beta(D, ativos, areasAtivas);

        if isempty(areasAtivas) && ~temRef
            % Sem lambda: Pg nominal
            tipo(ativos) = 2;
            ctrl = [];
            Cb   = sparse(N, 0);
        else
            usados = find(any(Cb ~= 0, 2));
            tipo(ativos) = 2;

            % GCIs com ajuste de P
            tipo(usados) = 3;
            ctrl.Cb      = Cb;
            ctrl.lig     = D.lig(areasAtivas);
            ctrl.NIesp   = D.NIesp(areasAtivas);

            % Areas com ramo no limite termico: a equacao NI_a = NIesp_a e
            % trocada pela equacao de fluxo ativo do ramo, com a mesma forma
            % de uma interligacao (r = terminal monitorado, t = outro).
            for j = 1:numel(areasAtivas)
                ia = areasAtivas(j);

                if termL(ia) > 0
                    ctrl.lig(j)  = CIA_LigRamo(D, termL(ia), termE(ia));
                    ctrl.NIesp(j) = -termP(ia);
                end
            end
            ctrl.temRef  = temRef;
            ctrl.PrefEsp = D.PrefEsp;
        end

        if isempty(Cb0)
            Cb0 = Cb;
        end


        %% d) Barras com Q Fixado no Limite
        idxQ = find(fixQ);

        for k = idxQ.'
            if tipo(k) == 3
                % GCI ativo com Q fixo
                tipo(k) = 4;
            else
                % PV -> PQ
                tipo(k) = 1;
            end
        end

        Q(idxQ) = QfixV(idxQ) - D.Qd(idxQ);


        %% e) Solucao do Estagio
        if strcmpi(cfg.partidaNR, 'quente')
            [V, lam, info] = CIA_Resolver(D, metodo, tipo, P, Q, ctrl, V0);
        else
            [V, lam, info] = CIA_Resolver(D, metodo, tipo, P, Q, ctrl);
        end

        if strcmpi(metodo, 'HELM')
            custo = custo + info.nrodadas;
        else
            custo = custo + info.n;
        end

        % Estagio sem solucao em gamma = 1
        if ~info.ok
            res = CIA_Montar(D, V, lam, info, tipo, fixP, PfixV, Cb, areasAtivas, loc, ativos);
            res.semSolucao = true;
            res.P = P;
            res.Q = Q;
            res.ctrl = ctrl;
            res.lamVet = lam;
            break
        end

        V0 = V;


        %% f) Retorno PQ -> PV (limites de Q)
        % Um gerador com Q no limite maximo cuja tensao ficou ACIMA do valor
        % especificado reduziria sua geracao de reativo para regular a tensao
        % (e o oposto no limite minimo). Nesses casos a barra volta a ter a
        % tensao controlada. cfg.maxRetornoQ limita os retornos de cada barra
        % para evitar oscilacao entre PV e PQ.
        if cfg.limQ && cfg.retornoQ
            kf  = find(fixQ);
            Vm  = abs(V(kf));
            noMax = QfixV(kf) >= D.Qmax(kf) - 1e-9;
            ret = (noMax & Vm > D.Vesp(kf) + cfg.tolVret) | ...
                  (~noMax & Vm < D.Vesp(kf) - cfg.tolVret);
            ret = ret & nRetQ(kf) < cfg.maxRetornoQ;

            if any(ret)
                kr = kf(ret);
                fixQ(kr)  = false;
                nRetQ(kr) = nRetQ(kr) + 1;
                regRet = [regRet; kr, est*ones(numel(kr),1)];
                continue
            end
        end


        %% g) Verificacao dos Limites (ordem definida por cfg.ordemLim)
        % 'QP': limites de Q primeiro, depois P dos GCIs
        % 'PQ': limites de P dos GCIs primeiro, depois Q
        if strcmpi(cfg.ordemLim, 'PQ')
            ordem = {'P', 'Q'};
        else
            ordem = {'Q', 'P'};
        end

        mudou = false;

        for o = 1:2

            % i. Limites de Potencia Reativa
            % cfg.nQporEstagio: numero maximo de barras convertidas por estagio
            % (Inf = todas as violadas; 1 = apenas a de maior violacao)
            if strcmp(ordem{o}, 'Q') && cfg.limQ
                [kq, ql, dq] = CIA_VerificaQ(D, V, tipo);

                if ~isempty(kq)
                    [~, ord] = sort(dq, 'descend');
                    nq  = min(numel(ord), cfg.nQporEstagio);
                    ord = ord(1:nq);
                    kq  = kq(ord);
                    ql  = ql(ord);
                    dq  = dq(ord);

                    fixQ(kq)  = true;
                    QfixV(kq) = ql;
                    regQ = [regQ; kq, ql, dq, est*ones(numel(kq),1)];
                    mudou = true;
                    break
                end
            end

            % ii. Limites de Potencia Ativa dos GCIs
            if strcmp(ordem{o}, 'P') && ~isempty(ctrl) && cfg.limP
                [kp, pl, dp] = CIA_VerificaP(D, ativos, Cb, lam);

                if ~isempty(kp)
                    if strcmpi(cfg.estrategia, 'um')
                        kp = kp(1);
                        pl = pl(1);
                        dp = dp(1);
                    end

                    fixP(kp)  = true;
                    PfixV(kp) = pl;
                    regP = [regP; kp, pl, dp, est*ones(numel(kp),1)];
                    mudou = true;
                    break
                end
            end
        end

        if mudou
            continue
        end


        %% h) Limite Termico dos Ramos
        if cfg.limTermico
            [termL, termE, termP, nTerm, regT, infoT, mudou, liberar] = ...
                CIA_VerificaTermico(D, V, lam, areasAtivas, termL, termE, ...
                                    termP, nTerm, regT, infoT, est);

            if liberar > 0
                fixP(D.area == liberar) = false;
            end

            if mudou
                continue
            end
        end


        %% i) Solucao Aceita
        res = CIA_Montar(D, V, lam, info, tipo, fixP, PfixV, Cb, areasAtivas, loc, ativos);
        res.semSolucao = false;
        res.P = P;
        res.Q = Q;
        res.ctrl = ctrl;
        res.lamVet = lam;
        break
    end



    %% ____________________________________________________________________
    %%
    %%                     3. INFORMACOES DO LACO
    %% ____________________________________________________________________

    res.regP      = regP;
    res.regQ      = regQ;
    res.Cb0       = Cb0;
    res.estagios  = est;
    res.reinicios = est - 1;
    res.custo     = custo;
    res.fixQ      = find(fixQ).';
    res.regRet    = regRet;
    res.termL     = termL;
    res.termE     = termE;
    res.termP     = termP;
    res.regT      = regT;
    res.infoT     = infoT;
end
