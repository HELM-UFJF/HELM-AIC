function [V, lam, info] = CIA_NR(Ybus, tipo, Vesp, P, Q, ref, ctrl, tol, maxit, V0)
%CIA_NR  Newton-Raphson com Controle de Intercambio de Area multiarea.
%
%   Resolve o mesmo problema de CIA_HELM, em coordenadas polares e com
%   Jacobiano esparso: as restricoes de intercambio NI_a = NIesp_a, medidas
%   nos terminais de cada area controlada, e um ajuste lambda_a por area
%   como incognita adicional. Se ctrl.temRef, a ultima coluna de Cb e a da
%   area da referencia, com a restricao P_ref = ctrl.PrefEsp.
%
%   Tipos de barra: 0 = referencia, 1 = PQ, 2 = PV,
%                   3 = GCI ativo (PV com P + dP), 4 = GCI ativo com Q no limite (PQ com P + dP),
%                   com dP_k = sum_a Cb(k,a)*lambda_a.
%   ctrl vazio => FP convencional.
%
%   V0 (opcional): tensoes complexas iniciais (partida quente); sem V0, partida plana.
%   Tudo em p.u.
%
%   Saidas: V, lambda (m x 1) e info (n: iteracoes, conv, erro).



    %% ____________________________________________________________________
    %%
    %%                       1. PREPARACAO
    %% ____________________________________________________________________

    %% a) Conjuntos de Barras e Posicoes no Vetor de Estado
    N = length(Vesp);

    if isempty(ctrl)
        m      = 0;
        mNI    = 0;
        temRef = false;
        Cb     = sparse(N, 0);
    else
        m      = size(ctrl.Cb, 2);

        % Restricoes de intercambio
        mNI    = numel(ctrl.lig);

        % Restricao de geracao da folga
        temRef = ctrl.temRef;
        Cb     = ctrl.Cb;
    end

    % Barras com equacao de Q
    pq  = find(tipo == 1 | tipo == 4);

    % Barras com equacao de P
    nr_ = find((1:N).' ~= ref);
    nth = numel(nr_);
    npq = numel(pq);
    nx  = nth + npq;

    posth = zeros(N,1);
    posth(nr_) = 1:nth;
    posv  = zeros(N,1);
    posv(pq)   = nth + (1:npq);


    %% b) Estimativa Inicial
    Vm  = Vesp;
    Vm(pq) = 1;
    th  = zeros(N,1);
    lam = zeros(m,1);

    if nargin >= 10 && ~isempty(V0)
        Vm(pq) = abs(V0(pq));
        th = angle(V0);
    end

    info.conv = false;
    info.erro = Inf;



    %% ____________________________________________________________________
    %%
    %%                        2. ITERACOES
    %% ____________________________________________________________________

    for it = 0:maxit

        %% a) Residuos
        V  = Vm.*exp(1j*th);
        I  = Ybus*V;
        Sc = V.*conj(I);

        Psp = P + Cb*lam;
        f   = [real(Sc(nr_)) - Psp(nr_); imag(Sc(pq)) - Q(pq)];

        for a = 1:mNI
            r  = ctrl.lig(a).r;
            t  = ctrl.lig(a).t;
            y  = ctrl.lig(a).y;
            g  = ctrl.lig(a).ysh;
            NI = -sum(real(V(r).*conj(y.*(V(r) - V(t)) + g.*V(r))));
            f  = [f; NI - ctrl.NIesp(a)];
        end

        if temRef
            f = [f; real(Sc(ref)) - ctrl.PrefEsp];
        end

        info.erro = max(abs(f));

        if info.erro < tol
            info.conv = true;
            break
        end


        %% b) Jacobiano Esparso
        dV   = spdiags(V, 0, N, N);
        dI   = spdiags(I, 0, N, N);
        dVn  = spdiags(V./abs(V), 0, N, N);
        dSth = 1j*dV*conj(dI - Ybus*dV);
        dSV  = dV*conj(Ybus*dVn) + conj(dI)*dVn;

        J = [real(dSth(nr_,nr_)), real(dSV(nr_,pq));
             imag(dSth(pq,nr_)),  imag(dSV(pq,pq))];


        %% c) Linhas e Colunas dos Intercambios
        if m > 0
            cols = sparse(nx, m);
            rows = zeros(m, nx);

            for a = 1:m
                ka = find(Cb(:,a) ~= 0);
                ka = ka(posth(ka) > 0);
                cols(:,a) = sparse(posth(ka), 1, -full(Cb(ka,a)), nx, 1);
            end

            if temRef
                rows(m,:) = [full(real(dSth(ref,nr_))), full(real(dSV(ref,pq)))];
            end

            for a = 1:mNI
                r = ctrl.lig(a).r;
                t = ctrl.lig(a).t;
                y = ctrl.lig(a).y;

                for i = 1:numel(r)
                    g  = real(y(i));
                    b  = imag(y(i));
                    d  = th(r(i)) - th(t(i));
                    c  = g*cos(d) + b*sin(d);
                    dc = -g*sin(d) + b*cos(d);

                    if posth(r(i)) > 0
                        rows(a,posth(r(i))) = rows(a,posth(r(i))) + Vm(r(i))*Vm(t(i))*dc;
                    end

                    if posth(t(i)) > 0
                        rows(a,posth(t(i))) = rows(a,posth(t(i))) - Vm(r(i))*Vm(t(i))*dc;
                    end

                    if posv(r(i)) > 0
                        rows(a,posv(r(i))) = rows(a,posv(r(i))) + (-2*g*Vm(r(i)) + Vm(t(i))*c);
                    end

                    if posv(t(i)) > 0
                        rows(a,posv(t(i))) = rows(a,posv(t(i))) + Vm(r(i))*c;
                    end

                    % Shunt do modelo pi: -real(ysh)*|V_r|^2
                    if posv(r(i)) > 0
                        rows(a,posv(r(i))) = rows(a,posv(r(i))) - 2*real(ctrl.lig(a).ysh(i))*Vm(r(i));
                    end
                end
            end

            J = [J, cols; sparse(rows), sparse(m, m)];
        end


        %% d) Atualizacao
        dx = J\(-f);
        th(nr_) = th(nr_) + dx(1:nth);
        Vm(pq)  = Vm(pq)  + dx(nth+1:nx);

        if m > 0
            lam = lam + dx(nx+1:end);
        end
    end

    info.n = it;
end
