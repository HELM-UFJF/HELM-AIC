function [V, lam, info] = CIA_Resolver(D, metodo, tipo, P, Q, ctrl, V0)
%CIA_RESOLVER  Resolve um estagio do controle pelo HELM ou pelo NR.
%
%   metodo: 'HELM' ou 'NR'. ctrl vazio => FP sem CIA.
%   V0 (opcional) e usado apenas pelo NR (partida quente).
%   info.ok indica se o estagio convergiu.


    %% a) Valor Padrao
    if nargin < 7
        V0 = [];
    end

    cfg = D.cfg;


    %% b) Chamada do Metodo
    if strcmpi(metodo, 'HELM')
        opc.xi     = cfg.xi;
        opc.nmax   = cfg.nmax;
        opc.nestag = cfg.nestag;
        opc.serie  = cfg.serie;
        opc.crit   = cfg.critHELM;

        [V, lam, info] = CIA_HELM(D.Ysr, D.Ysh, tipo, D.Vesp, P, Q, D.ref, ctrl, opc);
    else
        [V, lam, info] = CIA_NR(D.Ybus, tipo, D.Vesp, P, Q, D.ref, ctrl, ...
            cfg.tolNR, cfg.maxitNR, V0);
    end

    info.ok = info.conv;
end
