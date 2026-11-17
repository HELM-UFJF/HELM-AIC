function [kviol, Plim, dviol] = CIA_VerificaP(D, ativos, Cb, lam)
%CIA_VERIFICAP  Verifica os limites de geracao ativa dos GCIs ativos
%               com lambda(1) e devolve as violacoes em ordem decrescente.
%
%   Pg_k = Pg_nom_k + sum_a Cb(k,a)*lambda_a.
%   Saidas: GCIs violados, limite violado de cada um [p.u.] e violacao [p.u.].


    %% a) Geracao Resultante dos GCIs Ativos
    ativos = ativos(:);
    dP     = Cb*lam;
    pg     = D.Pg_nom(ativos) + dP(ativos);


    %% b) Violacoes
    acima  = pg - D.Pmax(ativos);
    abaixo = D.Pmin(ativos) - pg;
    dv     = max([acima, abaixo, zeros(numel(ativos),1)], [], 2);
    lim    = D.Pmin(ativos);
    lim(acima > 0) = D.Pmax(ativos(acima > 0));


    %% c) Ordenacao (maior violacao primeiro)
    sel = find(dv > D.cfg.epsP);
    [~, ord] = sort(dv(sel), 'descend');
    sel = sel(ord);

    kviol = ativos(sel);
    Plim  = lim(sel);
    dviol = dv(sel);
end
