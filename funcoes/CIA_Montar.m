function res = CIA_Montar(D, V, lam, info, tipo, fixP, PfixV, Cb, areasAtivas, loc, ativos)
%CIA_MONTAR  Organiza os resultados de um estagio aceito do CIA_Controle.
%
%   Inclui o despacho dos GCIs, os fatores de participacao realizados,
%   lambda de cada area controlada, a geracao da folga e os intercambios.


    %% a) Ajuste de Cada Area Controlada (NaN = perda de controlabilidade)
    lamA   = nan(D.nA, 1);
    lamRef = NaN;
    nAt    = numel(areasAtivas);

    if ~isempty(lam) && all(isfinite(lam))
        lamA(areasAtivas) = lam(1:nAt);
        dP = Cb*lam;

        if numel(lam) > nAt
            % Ajuste da area da folga ('folga')
            lamRef = lam(end);
        end
    else
        dP = zeros(D.N, 1);
    end


    %% b) Despacho dos GCIs [p.u.]
    Pg = D.Pg_nom;

    for k = D.gci(:).'
        if fixP(k)
            Pg(k) = PfixV(k);
        else
            Pg(k) = D.Pg_nom(k) + dP(k);
        end
    end


    %% c) Participacao Realizada na Variacao de Cada Area
    % Sinal positivo: a geracao do GCI aumentou; negativo: diminuiu.
    bf = nan(D.N, 1);

    for a = unique(D.area(D.gci)).'
        ks  = D.gci(D.area(D.gci) == a);
        tot = sum(Pg(ks) - D.Pg_nom(ks));

        if abs(tot) > 1e-14
            bf(ks) = (Pg(ks) - D.Pg_nom(ks))/abs(tot);
        end
    end


    %% d) Geracao da Folga e Reativo dos Geradores
    S      = V.*conj(D.Ybus*V);
    Pslack = real(S(D.ref)) + D.Pd(D.ref);
    Qg     = imag(S) + D.Qd;


    %% e) Estrutura de Saida
    res.V          = V;
    res.lam        = lamA;
    res.lamRef     = lamRef;
    res.info       = info;
    res.tipo       = tipo;
    res.loc        = loc;
    res.locAlguma  = any(loc);
    res.fixados    = find(fixP).';
    res.ativos     = ativos(:).';
    res.Cb         = Cb;
    res.beta_final = bf;
    res.Pg         = Pg;
    res.Qg         = Qg;
    res.Pslack     = Pslack;
    res.NI         = CIA_Intercambio(D, V);

    a1 = ativos(D.area(ativos) == D.areaRef);
    res.areaRefEsgotada = isempty(a1);
end
