function [V1, lam1, info] = CIA_HELM(Ysr, Ysh, tipo, Vesp, P, Q, ref, ctrl, opc)
%CIA_HELM  HELM com Controle de Intercambio de Area (CIA) multiarea embutido na recursao.
%
%   Resolve, para gamma = 1, o fluxo de potencia com as restricoes de
%   intercambio liquido NI_a = NIesp_a de cada area controlada embutidas, por
%   series de Maclaurin em gamma. A matriz de recursao, de dimensao 2N + m
%   (mais 2 com a regulacao da area de referencia), em que m e o numero de
%   areas controladas ativas, e esparsa e constante em n: e fatorada uma
%   unica vez e so o lado direito e atualizado a cada ordem.
%
%   Avaliacao das series em gamma = 1 (opc.serie):
%     'Wynn'        : algoritmo epsilon de Wynn incremental (Wynn_Passo); a
%                     cada ordem so a nova antidiagonal da tabela e calculada
%     'WynnCompleto': tabela epsilon refeita a cada ordem (Wynn_Epsilon),
%                     com as mesmas estimativas
%
%   Tipos de barra:
%     0 = referencia, 1 = PQ, 2 = PV,
%     3 = GCI ativo com tensao controlada (barra tipo H),
%     4 = GCI ativo com Q no limite (PQ com ajuste de potencia ativa)
%
%   Ajuste de potencia ativa dos GCIs: dP_k = sum_a Cb(k,a)*lambda_a (ver CIA_Beta).
%   Se ctrl.temRef, a ultima coluna de Cb e a da area da folga, com a
%   restricao P_ref = ctrl.PrefEsp embutida (equacao real da barra de
%   referencia), e Q_ref passa a ser uma incognita da recursao.
%
%   Entradas (tudo em p.u.):
%     Ysr, Ysh : matriz nodal sem shunts (esparsa) e vetor de shunts
%     Vesp     : modulo de tensao especificado (PV, H e referencia)
%     P, Q     : injecoes liquidas especificadas (nas barras 3 e 4, P e o valor nominal)
%     ctrl     : estrutura do controle (vazia => FP sem CIA)
%                ctrl.Cb    : N x m, coeficientes de participacao
%                ctrl.lig(a): linhas de intercambio da area a (campos r, t, y)
%                ctrl.NIesp : m x 1, intercambio liquido recebido especificado
%     opc      : estrutura com xi, nmax, nestag e serie
%
%   Saidas: tensoes complexas V1, lambda(1) de cada area (m x 1) e info
%           (n: ordem final, conv, erro, hist, estag).



    %% ____________________________________________________________________
    %%
    %%                       1. PREPARACAO
    %% ____________________________________________________________________

    %% a) Conjuntos de Barras
    N     = length(Ysh);

    % Barras com |V| especificado
    PVH   = find(tipo == 2 | tipo == 3);
    isPVH = false(N,1);
    isPVH(PVH) = true;

    if isempty(ctrl)
        m = 0;
    else
        m = size(ctrl.Cb, 2);
    end

    if m > 0
        Cb  = ctrl.Cb;

        % Barras com ajuste de P
        GCI = find(any(Cb ~= 0, 2));
    else
        Cb  = sparse(N, 0);
        GCI = [];
    end

    if m > 0
        % Restricoes de intercambio
        mNI    = numel(ctrl.lig);

        % Restricao de geracao da folga
        temRef = ctrl.temRef;
    else
        mNI    = 0;
        temRef = false;
    end

    % + Q_ref quando ha restricao da folga
    nu   = 2*N + m + temRef;

    % Menor ordem com avaliacao em gamma = 1
    nmin = 4;


    %% b) Elementos da Matriz de Recursao (formato de coordenadas)
    % Linhas 2k-1 (parte real) e 2k (parte imaginaria) de sum_j Ysr(k,j) V_j[n]
    [ii, jj, yy] = find(Ysr);
    gg = real(yy);
    bb = imag(yy);

    Ia = [2*ii-1; 2*ii-1; 2*ii; 2*ii];
    Ja = [2*jj-1; 2*jj;   2*jj-1; 2*jj];
    Va = [gg;     -bb;    bb;     gg];


    %% c) Barras PV/H: a coluna de Vre vira a coluna de Q
    colPV = false(nu, 1);
    colPV(2*PVH-1) = true;
    manter = ~colPV(Ja);
    Ia = [Ia(manter); 2*PVH];
    Ja = [Ja(manter); 2*PVH-1];
    Va = [Va(manter); ones(numel(PVH),1)];


    %% d) Colunas de lambda (equacao real das barras GCI)
    for a = 1:m
        ka = find(Cb(:,a) ~= 0);
        Ia = [Ia; 2*ka-1];
        Ja = [Ja; (2*N+a)*ones(numel(ka),1)];
        Va = [Va; -full(Cb(ka,a))];
    end


    %% e) Barra de Referencia: V fixo
    manter = (Ia ~= 2*ref-1) & (Ia ~= 2*ref);
    Ia = [Ia(manter); 2*ref-1; 2*ref];
    Ja = [Ja(manter); 2*ref-1; 2*ref];
    Va = [Va(manter); 1; 1];


    %% f) Linhas das Restricoes de Intercambio (uma por area)
    cre = zeros(N, mNI);
    cim = zeros(N, mNI);

    for a = 1:mNI
        r = ctrl.lig(a).r;
        t = ctrl.lig(a).t;
        y = ctrl.lig(a).y;

        for i = 1:numel(r)
            cre(r(i),a) = cre(r(i),a) + real(y(i));
            cre(t(i),a) = cre(t(i),a) - real(y(i));
            cim(r(i),a) = cim(r(i),a) - imag(y(i));
            cim(t(i),a) = cim(t(i),a) + imag(y(i));
        end

        kre = find(cre(:,a) ~= 0 & ~isPVH);
        kim = find(cim(:,a) ~= 0);

        Ia = [Ia; (2*N+a)*ones(numel(kre),1); (2*N+a)*ones(numel(kim),1)];
        Ja = [Ja; 2*kre-1; 2*kim];
        Va = [Va; cre(kre,a); cim(kim,a)];
    end


    %% g) Restricao da Area da Folga: P_ref = PrefEsp
    % Linha 2N+m : parte real da equacao da barra de referencia (coluna de
    %              lambda_ref ja incluida pela ultima coluna de Cb)
    % Linha 2N+m+1: parte imaginaria, com a incognita Q_ref (coluna 2N+m+1)
    if temRef
        [~, jr, yr] = find(Ysr(ref,:));
        jr  = jr(:);
        yr  = yr(:);
        gr  = real(yr);
        br  = imag(yr);
        nPV = ~isPVH(jr);
        lR  = 2*N + m;
        lI  = 2*N + m + 1;

        Ia = [Ia; lR*ones(sum(nPV),1); lR*ones(numel(jr),1); ...
                  lI*ones(sum(nPV),1); lI*ones(numel(jr),1); lI];
        Ja = [Ja; 2*jr(nPV)-1; 2*jr; 2*jr(nPV)-1; 2*jr; lI];
        Va = [Va; gr(nPV); -br; br(nPV); gr; 1];
    end


    %% h) Fatoracao LU Esparsa (unica)
    % Forma com vetores de permutacao: A(p,q) = L*U
    A = sparse(Ia, Ja, Va, nu, nu);
    [L, U, pL, qL] = lu(A, 'vector');


    %% i) Termos Constantes Usados em Todas as Ordens
    % Evita extrair submatrizes esparsas dentro do laco das ordens.
    % Acoplamento de Vre das barras PV/H
    GPV  = real(Ysr(:,PVH));
    BPV  = imag(Ysr(:,PVH));

    % Coeficientes dos GCIs (denso)
    CbG  = full(Cb(GCI,:));

    % Parcela de Vre nas restricoes NI
    creP = cre(PVH,:);

    if temRef
        GPVr = full(real(Ysr(ref,PVH)));
        BPVr = full(imag(Ysr(ref,PVH)));

        % Linha da referencia (extraida uma vez)
        Yref = Ysr(ref,:);
    end



    %% ____________________________________________________________________
    %%
    %%                     2. RECURSAO DAS SERIES
    %% ____________________________________________________________________

    %% a) Solucao Germinal (coluna n+1 guarda a ordem n)
    nmax = opc.nmax;
    V    = zeros(N, nmax+1);
    W    = zeros(N, nmax+1);
    Qh   = zeros(N, nmax+1);
    lam  = zeros(m, nmax+1);
    QI   = zeros(mNI, nmax+1);

    % Serie de Q da barra de referencia
    qref = zeros(1, nmax+1);

    V(:,1) = 1;
    W(:,1) = 1;
    Sesp   = P - 1j*Q;


    %% b) Controle de Convergencia
    hist    = [];
    conv    = false;
    estag   = false;
    Vp_old  = [];
    QIp_old = zeros(mNI,1);
    Vp      = V(:,1);
    lamp    = zeros(m,1);
    errbest = Inf;
    nbest   = 0;
    Vbest   = Vp;
    lambest = zeros(m,1);
    n       = 0;


    % Criterio de parada
    resid = isfield(opc, 'crit') && strcmpi(opc.crit, 'residuo');

    if resid
        Ybus = Ysr + spdiags(Ysh, 0, N, N);
        nRef = setdiff((1:N).', ref);
        isPQ = (tipo == 1 | tipo == 4);
    end

    % Wynn incremental: estados das tabelas epsilon (soma parcial de ordem 0)
    incr = strcmpi(opc.serie, 'Wynn');

    % As series de V, QI e lambda sao empilhadas numa unica tabela (o
    % algoritmo e independente linha a linha): uma chamada por ordem.
    if incr
        iV = 1:N;
        iQ = N + (1:mNI);
        iL = N + mNI + (1:m);
        sS = [V(:,1); QI(:,1); lam(:,1)];
        [wS, ~] = Wynn_Passo([], sS);
    end


    %% c) Laco das Ordens
    for n = 1:nmax

        % i. Parte real conhecida das barras PV/H
        Vre = zeros(N,1);

        if n == 1
            Vre(PVH) = (Vesp(PVH).^2 - 1)/2;
        else
            Vre(PVH) = -0.5*real(sum(V(PVH,2:n).*conj(V(PVH,n:-1:2)), 2));
        end

        % ii. Lado direito
        R = Sesp.*conj(W(:,n)) - Ysh.*V(:,n);

        if n > 1
            % So barras PV/H
            sq = sum(Qh(PVH,2:n).*conj(W(PVH,n:-1:2)), 2);
        else
            sq = zeros(numel(PVH),1);
        end

        R(PVH) = P(PVH).*conj(W(PVH,n)) - 1j*sq - Ysh(PVH).*V(PVH,n);

        if m > 0 && n > 1
            % sum_a Cb(k,a) * sum_{j=1}^{n-1} lambda_a[j] W_k*[n-j]
            % |GCI| x m
            conv_lw = conj(W(GCI,n:-1:2))*lam(:,2:n).';
            R(GCI)  = R(GCI) + sum(CbG.*conv_lw, 2);
        end

        rhs = zeros(nu,1);
        rhs(1:2:2*N) = real(R) - GPV*Vre(PVH);
        rhs(2:2:2*N) = imag(R) - BPV*Vre(PVH);
        rhs(2*ref-1) = (n == 1)*(Vesp(ref) - 1);
        rhs(2*ref)   = 0;

        for a = 1:mNI
            r = ctrl.lig(a).r;
            t = ctrl.lig(a).t;
            y = ctrl.lig(a).y;
            g = ctrl.lig(a).ysh;
            s = 0;

            if n > 1
                s = real(sum(conj(y).*sum(V(r,2:n).*(conj(V(r,n:-1:2)) - conj(V(t,n:-1:2))), 2)));
            end

            % Shunt do modelo pi (embutido com gamma): ordem n-1 de V_r V_r*
            sh = real(sum(conj(g).*sum(V(r,1:n).*conj(V(r,n:-1:1)), 2)));

            rhs(2*N+a) = -(n == 1)*ctrl.NIesp(a) - s - sh - sum(creP(:,a).*Vre(PVH));
        end

        if temRef
            Rr = ctrl.PrefEsp*conj(W(ref,n)) - Ysh(ref)*V(ref,n);

            if n > 1
                Rr = Rr - 1j*sum(qref(2:n).*conj(W(ref,n:-1:2)));
            end

            rhs(2*N+m)   = real(Rr) - GPVr*Vre(PVH);
            rhs(2*N+m+1) = imag(Rr) - BPVr*Vre(PVH);
        end

        % iii. Solucao do sistema linear da ordem n
        x = zeros(nu,1);
        x(qL) = U\(L\rhs(pL));
        Vn = x(1:2:2*N) + 1j*x(2:2:2*N);
        Vn(PVH) = Vre(PVH) + 1j*x(2*PVH);

        V(:,n+1)    = Vn;
        Qh(PVH,n+1) = x(2*PVH-1);

        if m > 0
            lam(:,n+1) = x(2*N+1:2*N+m);
        end

        if temRef
            qref(n+1) = x(2*N+m+1);
        end

        % iv. Reciproco da tensao (convolucao): W[n] = -sum_{j=0}^{n-1} W[j] V[n-j]
        W(:,n+1) = -sum(W(:,1:n).*V(:,n+1:-1:2), 2);

        % v. Potencia reativa de intercambio da ordem n (por area)
        for a = 1:mNI
            r = ctrl.lig(a).r;
            t = ctrl.lig(a).t;
            y = ctrl.lig(a).y;
            g = ctrl.lig(a).ysh;
            QI(a,n+1) = -imag(sum(conj(y).*sum(V(r,1:n+1).*(conj(V(r,n+1:-1:1)) - conj(V(t,n+1:-1:1))), 2))) ...
                        - imag(sum(conj(g).*sum(V(r,1:n).*conj(V(r,n:-1:1)), 2)));
        end

        % vi. Avaliacao em gamma = 1
        % 'Wynn': forma incremental (so a nova antidiagonal a cada ordem);
        % 'WynnCompleto': tabela refeita a cada ordem (mesmo resultado)
        if incr
            sS = sS + [V(:,n+1); QI(:,n+1); lam(:,n+1)];
            [wS, eS] = Wynn_Passo(wS, sS);
            Vp   = eS(iV);
            QIp  = real(eS(iQ));
            lamp = real(eS(iL));
        end

        if n < nmin
            hist(end+1) = Inf;
            continue
        end

        if ~incr
            Vp   = Wynn_Epsilon(cumsum(V(:,1:n+1), 2));
            QIp  = real(Wynn_Epsilon(cumsum(QI(:,1:n+1), 2)));
            lamp = real(Wynn_Epsilon(cumsum(lam(:,1:n+1), 2)));
        end

        if isempty(Vp_old)
            Vp_old  = Vp;
            QIp_old = QIp;
            hist(end+1) = Inf;
            continue
        end

        % vii. Criterio de convergencia (tensoes, intercambios e geracao da referencia)
        % 'passo'  : diferenca entre estimativas sucessivas de V e QI;
        % 'residuo': maior desbalanco de P, Q e |V| na estimativa atual, o
        %            mesmo criterio do NR (cfg.tolNR = cfg.xi)
        if resid
            Sc = Vp.*conj(Ybus*Vp);
            Pe = P;

            if m > 0
                Pe = Pe + Cb*lamp;
            end

            dP = abs(real(Sc(nRef)) - Pe(nRef));
            dQ = abs(imag(Sc(isPQ)) - Q(isPQ));
            dVm = abs(abs(Vp(PVH)) - Vesp(PVH));
            e1 = max([0; dP; dQ; dVm]);
        else
            e1 = max(abs(Vp - Vp_old));
        end

        e2 = 0;

        for a = 1:mNI
            r  = ctrl.lig(a).r;
            t  = ctrl.lig(a).t;
            y  = ctrl.lig(a).y;
            g  = ctrl.lig(a).ysh;
            NI = -sum(real(Vp(r).*conj(y.*(Vp(r) - Vp(t)) + g.*Vp(r))));
            e2 = max(e2, abs(NI - ctrl.NIesp(a)));
        end

        if temRef
            Pr = real(Vp(ref)*conj(Yref*Vp + Ysh(ref)*Vp(ref)));
            e2 = max(e2, abs(Pr - ctrl.PrefEsp));
        end

        if resid
            % NI ja verificado em e2
            e3 = 0;
        else
            e3 = max([0; abs(QIp - QIp_old)]);
        end

        err = max([e1, e2, e3]);
        hist(end+1) = err;

        Vp_old  = Vp;
        QIp_old = QIp;

        if err < opc.xi
            conv = true;
            break
        end

        % viii. Parada por estagnacao
        if err < errbest
            errbest = err;
            nbest   = n;
            Vbest   = Vp;
            lambest = lamp;
        elseif n - nbest >= opc.nestag
            estag = true;
            break
        end
    end



    %% ____________________________________________________________________
    %%
    %%                          3. SAIDAS
    %% ____________________________________________________________________

    %% a) Solucao
    if conv || nbest == 0
        V1   = Vp;
        lam1 = real(lamp);
        nfim = n;
    else
        % Nao convergiu: ordem de menor erro
        V1   = Vbest;
        lam1 = real(lambest);
        nfim = nbest;
    end


    %% b) Informacoes
    info.n        = nfim;
    info.nrodadas = n;
    info.conv     = conv;
    info.estag    = estag;
    info.erro     = min([errbest, hist(end)]);
    info.hist     = hist;
end
