function [kviol, Qlim, dviol] = CIA_VerificaQ(D, V, tipo)
%CIA_VERIFICAQ  Verifica os limites de potencia reativa das barras com tensao
%               controlada (PV e H), exceto a referencia.
%
%   Qg_k = Q_liquido_k + Qd_k deve satisfazer Qmin_k <= Qg_k <= Qmax_k.
%   Saidas: barras violadas (todas), limite violado de cada uma [p.u.] e
%           o tamanho da violacao [p.u.].


    %% a) Potencia Reativa Gerada
    S  = V.*conj(D.Ybus*V);
    Qg = imag(S) + D.Qd;


    %% b) Candidatas: PV e H
    cand = find(tipo == 2 | tipo == 3);


    %% c) Violacoes
    acima  = Qg(cand) - D.Qmax(cand);
    abaixo = D.Qmin(cand) - Qg(cand);
    dv     = max([acima, abaixo, zeros(numel(cand),1)], [], 2);
    sel    = dv > D.cfg.epsQ;

    kviol = cand(sel);
    dviol = dv(sel);
    Qlim  = zeros(numel(kviol), 1);

    for i = 1:numel(kviol)
        k = kviol(i);

        if Qg(k) > D.Qmax(k)
            Qlim(i) = D.Qmax(k);
        else
            Qlim(i) = D.Qmin(k);
        end
    end
end
