function [dmax, d] = CIA_Residuo(D, res)
%CIA_RESIDUO  Maior residuo (mismatch) das equacoes do ultimo estagio [p.u.].
%
%   Avalia, na solucao (V, lambda) devolvida por CIA_Controle, os residuos:
%     P  : todas as barras exceto a referencia (com o ajuste Cb*lambda)
%     Q  : barras PQ (tipos 1 e 4)
%     NI : restricoes de intercambio das areas controladas ativas
%     Pref: geracao da barra de referencia (modo 'folga')
%   E a mesma medida para o HELM e para o NR, como na Fig. 6 de Luo e Liu.
%
%   Saidas: dmax = maior residuo em modulo; d = estrutura com cada parcela.


    %% a) Injecoes Calculadas
    V  = res.V;
    S  = V.*conj(D.Ybus*V);
    P  = res.P;
    Q  = res.Q;
    lam = res.lamVet;
    ctrl = res.ctrl;

    if ~isempty(ctrl) && ~isempty(lam) && all(isfinite(lam))
        P = P + ctrl.Cb*lam;
    end


    %% b) Residuos de P e Q
    nr = (1:D.N).' ~= D.ref;
    pq = (res.tipo == 1 | res.tipo == 4);
    d.P = max(abs(real(S(nr)) - P(nr)));
    d.Q = max([0; abs(imag(S(pq)) - Q(pq))]);


    %% c) Residuos das Restricoes de Intercambio e da Folga
    d.NI  = 0;
    d.Ref = 0;

    if ~isempty(ctrl)
        for a = 1:numel(ctrl.lig)
            r  = ctrl.lig(a).r;
            t  = ctrl.lig(a).t;
            y  = ctrl.lig(a).y;
            g  = ctrl.lig(a).ysh;
            NI = -sum(real(V(r).*conj(y.*(V(r) - V(t)) + g.*V(r))));
            d.NI = max(d.NI, abs(NI - ctrl.NIesp(a)));
        end

        if ctrl.temRef
            d.Ref = abs(real(S(D.ref)) - ctrl.PrefEsp);
        end
    end


    %% d) Maior Residuo
    dmax = max([d.P, d.Q, d.NI, d.Ref]);

    if ~all(isfinite(V))
        dmax = Inf;
    end
end
