function lig = CIA_LigRamo(D, L, e)
%CIA_LIGRAMO  Estrutura "lig" de um unico ramo, com terminal monitorado e.
%
%   O fluxo ativo que SAI da barra e para o ramo L e P = Re{V_e conj(I_e)}.
%   Com lig.r = e, a funcao de intercambio retorna NI = -P, de modo que a
%   restricao NI = -P_alvo impoe P = P_alvo no terminal monitorado.


    %% a) Terminais e Admitancias do Modelo pi
    r = D.ramo;

    if e == r.de(L)
        lig.r   = r.de(L);
        lig.t   = r.para(L);
        lig.ysh = r.yshDe(L);
    else
        lig.r   = r.para(L);
        lig.t   = r.de(L);
        lig.ysh = r.yshPara(L);
    end

    lig.y = r.y(L);
end
