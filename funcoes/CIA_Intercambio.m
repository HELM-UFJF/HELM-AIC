function NI = CIA_Intercambio(D, V)
%CIA_INTERCAMBIO  Intercambio liquido recebido por cada area controlada [p.u.].
%
%   NI(ia) = soma dos fluxos de potencia ativa nas linhas de intercambio da
%   area D.areasCtrl(ia), medidos nos terminais r dentro da area,
%   com o modelo pi completo do ramo (inclusive o tap dos transformadores).
%   Negativo: a area exporta.


    %% a) Fluxo Liquido por Area
    NI = zeros(D.nA, 1);

    for ia = 1:D.nA
        r = D.lig(ia).r;
        t = D.lig(ia).t;
        y = D.lig(ia).y;
        g = D.lig(ia).ysh;
        NI(ia) = -sum(real(V(r).*conj(y.*(V(r) - V(t)) + g.*V(r))));
    end
end
