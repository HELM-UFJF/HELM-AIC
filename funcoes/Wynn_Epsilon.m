function est = Wynn_Epsilon(S)
%WYNN_EPSILON  Algoritmo epsilon de Wynn aplicado as somas parciais de varias
%              series ao mesmo tempo (uma serie por linha de S).
%
%   S(:,k+1) = soma parcial de ordem k (k = 0, ..., n) de cada serie.
%   As colunas pares da tabela epsilon coincidem com aproximantes de Pade
%   (diagonais e quase diagonais) avaliados em gamma = 1, mas sao obtidas
%   apenas com somas e divisoes elemento a elemento, sem resolver sistemas
%   lineares nem calcular determinantes. Por isso o calculo e vetorizado
%   para todas as barras de uma vez.
%
%   Recursao:  eps_{-1}^(k) = 0,  eps_0^(k) = S_k,
%              eps_{j+1}^(k) = eps_{j-1}^(k+1) + 1/(eps_j^(k+1) - eps_j^(k))
%
%   Saida: estimativa do limite de cada serie (ultimo elemento da coluna par
%          mais alta da tabela).


    %% a) Inicializacao da Tabela
    [nl, m] = size(S);

    % Coluna eps_{-1}
    eAnt = zeros(nl, m+1);

    % Coluna eps_0
    eAtu = S;
    est  = S(:,end);

    % Linhas em que a tabela segue valida
    ativa = true(nl, 1);
    j    = 0;


    %% b) Construcao das Colunas
    while size(eAtu, 2) > 1 && any(ativa)
        L = size(eAtu, 2);
        d = eAtu(:,2:L) - eAtu(:,1:L-1);

        % i. Quebra da tabela: diferenca desprezivel (serie ja estacionaria
        %    ou divisao por zero). A linha para e mantem a ultima estimativa.
        escala  = max(abs(eAtu(:,2:L)), abs(eAtu(:,1:L-1)));
        pequeno = abs(d) <= 1e-13*max(escala, 1e-300);
        ativa   = ativa & ~any(pequeno, 2);

        % Valor inocuo (linha ja inativa)
        d(pequeno) = 1;

        % ii. Nova coluna
        eNov = eAnt(:,2:L) + 1./d;
        eAnt = eAtu;
        eAtu = eNov;
        j    = j + 1;

        % iii. Colunas pares aproximam o limite
        if mod(j, 2) == 0
            cand = eAtu(:,end);
            ok   = ativa & isfinite(cand);
            est(ok) = cand(ok);
        end
    end
end
