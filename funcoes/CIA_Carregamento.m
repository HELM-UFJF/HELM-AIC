function [carg, Smax] = CIA_Carregamento(mpc, V, rateA)
%CIA_CARREGAMENTO  Carregamento dos ramos [p.u. da capacidade] (modelo pi com tap).
%
%   carg(L) = max(|S_de|, |S_para|) / rateA(L), com S em MVA.


    %% a) Correntes nos Terminais
    br = mpc.branch;
    f  = br(:,1);
    t  = br(:,2);
    a  = br(:,9);
    a(a == 0) = 1;
    y  = 1./(br(:,3) + 1j*br(:,4));
    bc = 1j*br(:,5)/2;

    If = (y./a.^2 + bc).*V(f) - (y./a).*V(t);
    It = (y + bc).*V(t) - (y./a).*V(f);


    %% b) Potencias e Carregamento
    Smax = mpc.baseMVA*max(abs(V(f).*conj(If)), abs(V(t).*conj(It)));
    carg = Smax./rateA;
end
