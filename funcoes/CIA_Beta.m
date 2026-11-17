function [Cb, temRef] = CIA_Beta(D, ativos, areasAtivas)
%CIA_BETA  Matriz de coeficientes de participacao Cb dos GCIs ativos.
%
%   Cada area controlada ativa a (coluna de Cb) tem um ajuste lambda_a e
%   o GCI k ativo da area a recebe Cb(k,a) = +beta_k, com
%   beta_k = Pg_nom_k / sum(Pg_nom dos GCIs ativos da mesma area).
%   Assim dP_k = sum_a Cb(k,a)*lambda_a.
%
%   Papel dos GCIs da area da folga (D.cfg.modoArea1):
%     'espelho': Cb(k,a) = -beta_k em todas as colunas, isto e, compensam
%                o ajuste de todas as areas (formulacao de duas areas,
%                em que lambda_2 = -lambda).
%     'folga'  : a area da folga tem o seu proprio ajuste lambda_ref (ultima
%                coluna de Cb), usado para manter a geracao da barra de
%                referencia no valor programado. E o equivalente ao CAG da
%                area da folga: ela cobre o proprio balanco com os seus GCIs,
%                e a barra de referencia fica apenas com o seu despacho.
%     'nenhum' : os GCIs da area da folga nao participam.
%
%   areasAtivas: indices (em D.areasCtrl) das areas controladas que ainda tem
%                GCI ativo.
%   temRef     : true se a ultima coluna de Cb e a da area da folga ('folga').


    %% a) Fatores Normalizados por Area
    beta = zeros(D.N, 1);

    for a = unique(D.area(ativos)).'
        act = ativos(D.area(ativos) == a);

        if D.cfg.betaIgual
            % Fatores iguais (Luo e Liu)
            beta(act) = 1/numel(act);
        else
            beta(act) = D.Pg_nom(act)/sum(D.Pg_nom(act));
        end
    end


    %% b) Colunas das Areas Controladas
    m  = numel(areasAtivas);
    a1 = ativos(D.area(ativos) == D.areaRef);
    temRef = strcmpi(D.cfg.modoArea1, 'folga') && ~isempty(a1);
    Cb = sparse(D.N, m + temRef);

    for j = 1:m
        a   = D.areasCtrl(areasAtivas(j));
        act = ativos(D.area(ativos) == a);
        Cb(act, j) = beta(act);

        if strcmpi(D.cfg.modoArea1, 'espelho')
            Cb(a1, j) = -beta(a1);
        end
    end


    %% c) Coluna da Area da Folga
    if temRef
        Cb(a1, m+1) = beta(a1);
    end
end
