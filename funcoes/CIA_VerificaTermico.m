function [termL, termE, termP, nTerm, regT, infoT, mudou, liberar] = ...
    CIA_VerificaTermico(D, V, lam, areasAtivas, termL, termE, termP, nTerm, regT, infoT, est)
%CIA_VERIFICATERMICO  Limite termico dos ramos no CIA (troca de equacao).
%
%   Apos um estagio convergido, verifica o carregamento de todos os ramos.
%   i.  Area ja limitada: o alvo de P do ramo monitorado e corrigido com o Q
%       atual, P_lim = sign(P)*sqrt(S_max^2 - Q^2), ate |S| = S_max.
%       Se outro ramo da mesma area ficar mais carregado, troca-se o ramo.
%   ii. Ramo violado sem controle: o ramo e atribuido a uma area controlada
%       ativa e a equacao NI dessa area e trocada pela de fluxo do ramo:
%         - ramo interno a uma area controlada -> essa area;
%         - interligacao com a area de referencia -> a area controlada;
%         - interligacao entre areas controladas -> a de maior |lambda|;
%         - ramo interno a area de referencia -> apenas aviso.
%   Uma alteracao por estagio (a de maior violacao).
%
%   liberar: area (codigo) cujos GCIs fixados no limite de P devem ser
%   liberados. Quando o ramo passa a definir o controle da area, o
%   intercambio programado deixa de valer e os GCIs que saturaram para
%   atende-lo voltam a participar (inclusive em area sem controle, LOC).


    %% a) Carregamento Atual
    cfg   = D.cfg;
    tol   = cfg.tolTerm;
    r     = D.ramo;
    Smax  = r.rate/D.baseMVA;
    [Sde, Spara, carga] = CIA_FluxoRamos(D, V);
    mudou   = false;
    liberar = 0;
    aDe   = D.area(r.de);
    aPara = D.area(r.para);


    %% b) Refinamento das Areas ja Limitadas
    for ia = find(termL > 0)
        L = termL(ia);

        if ~any(areasAtivas == ia)
            msg = sprintf('Area %d sem GCI livre: ramo %d-%d nao controlado', ...
                          D.areasCtrl(ia), r.de(L), r.para(L));

            if ~any(strcmp(infoT, msg))
                infoT{end+1} = msg;
            end

            % Area sem controle (LOC)
            continue
        end

        % Ramo mais carregado da area (exceto ja monitorado) -> troca
        naArea = find((aDe == D.areasCtrl(ia) | aPara == D.areasCtrl(ia)) & carga > 1 + tol);
        naArea = setdiff(naArea, L);
        naArea = setdiff(naArea, termL(termL > 0 & (1:D.nA) ~= ia));

        if ~isempty(naArea) && max(carga(naArea)) > carga(L) && nTerm(ia) < cfg.maxTerm
            [~, k] = max(carga(naArea));
            L2 = naArea(k);
            [termE(ia), termP(ia)] = alvo(r, L2, Sde, Spara, Smax);
            termL(ia) = L2;
            nTerm(ia) = nTerm(ia) + 1;
            liberar   = D.areasCtrl(ia);
            regT = [regT; D.areasCtrl(ia), L2, termP(ia)*D.baseMVA, est];
            mudou = true;
            return
        end

        % Correcao do alvo de P com o Q atual
        if abs(carga(L) - 1) > tol && nTerm(ia) < cfg.maxTerm
            [termE(ia), termP(ia), ok] = alvo(r, L, Sde, Spara, Smax);

            if ~ok
                infoT{end+1} = sprintf('Ramo %d-%d: Q excede a capacidade', r.de(L), r.para(L));
            end

            nTerm(ia) = nTerm(ia) + 1;
            regT = [regT; D.areasCtrl(ia), L, termP(ia)*D.baseMVA, est];
            mudou = true;
            return
        end
    end


    %% c) Novos Ramos Violados
    viol = find(carga > 1 + tol);
    viol = setdiff(viol, termL(termL > 0));
    [~, ord] = sort(carga(viol), 'descend');
    viol = viol(ord);

    for k = 1:numel(viol)
        L  = viol(k);
        ia = areaDoRamo(D, aDe(L), aPara(L), lam, areasAtivas);

        if ia == 0
            msg = sprintf('Ramo %d-%d (%.1f %%) sem area controlada ativa', ...
                          r.de(L), r.para(L), 100*carga(L));

            if ~any(strcmp(infoT, msg))
                infoT{end+1} = msg;
            end

            continue
        end

        if termL(ia) > 0 || nTerm(ia) >= cfg.maxTerm
            % Ja tratado em b)
            continue
        end

        [termE(ia), termP(ia)] = alvo(r, L, Sde, Spara, Smax);
        termL(ia) = L;
        nTerm(ia) = nTerm(ia) + 1;
        liberar   = D.areasCtrl(ia);
        regT = [regT; D.areasCtrl(ia), L, termP(ia)*D.baseMVA, est];
        mudou = true;
        return
    end
end



%% ________________________________________________________________________
%%
%%                        FUNCOES AUXILIARES
%% ________________________________________________________________________

function [e, Pa, ok] = alvo(r, L, Sde, Spara, Smax)
%ALVO  Terminal de maior |S| e fluxo ativo alvo com o Q atual.

    if abs(Sde(L)) >= abs(Spara(L))
        e = r.de(L);
        S = Sde(L);
    else
        e = r.para(L);
        S = Spara(L);
    end

    q2 = Smax(L)^2 - imag(S)^2;
    ok = q2 > 0;
    Pa = sign(real(S))*sqrt(max(q2, 0));
end


function ia = areaDoRamo(D, a1, a2, lam, areasAtivas)
%AREADORAMO  Area controlada (indice em areasCtrl) que assume o ramo.
%   Candidatas: areas controladas com GCIs (livres ou fixados). Entre duas
%   candidatas ativas, a de maior |lambda|; se uma estiver sem controle
%   (GCIs saturados), ela e a que empurra o fluxo e assume o ramo.

    ia = 0;
    i1 = find(D.areasCtrl == a1);
    i2 = find(D.areasCtrl == a2);
    i1 = i1(any(D.area(D.gci) == a1));
    i2 = i2(any(D.area(D.gci) == a2));

    if isempty(i1) && isempty(i2)
        return
    elseif isempty(i2)
        ia = i1;
    elseif isempty(i1)
        ia = i2;
    elseif i1 == i2
        ia = i1;
    else
        j1 = find(areasAtivas == i1);
        j2 = find(areasAtivas == i2);

        if isempty(j1)
            ia = i1;
        elseif isempty(j2)
            ia = i2;
        elseif abs(lam(j1)) >= abs(lam(j2))
            ia = i1;
        else
            ia = i2;
        end
    end
end
