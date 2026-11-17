function [Sde, Spara, carga] = CIA_FluxoRamos(D, V)
%CIA_FLUXORAMOS  Fluxo de potencia complexa nos ramos e carregamento.
%
%   Sde, Spara : potencia que SAI de cada terminal para o ramo [p.u.]
%   carga      : max(|Sde|, |Spara|)*baseMVA / capacidade (rateA)


    %% a) Correntes nos Terminais (modelo pi com tap)
    r  = D.ramo;
    Vf = V(r.de);
    Vt = V(r.para);
    If = r.y.*(Vf - Vt) + r.yshDe.*Vf;
    It = r.y.*(Vt - Vf) + r.yshPara.*Vt;


    %% b) Potencias e Carregamento
    Sde   = Vf.*conj(If);
    Spara = Vt.*conj(It);
    carga = D.baseMVA*max(abs(Sde), abs(Spara))./r.rate;
end
