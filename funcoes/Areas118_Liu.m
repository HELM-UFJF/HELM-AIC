function area = Areas118_Liu()
%AREAS118_LIU  Nove areas do IEEE 118 barras de Luo e Liu (EPSR 2023, Fig. 4).
%
%   A divisao nao e listada no artigo; foi reconstruida a partir da Fig. 4 e
%   conferida pelos intercambios originais da Tabela 3 (case118 do MATPOWER,
%   fluxos de linhas e transformadores, exportacao positiva):
%     Areas 1, 2, 4, 5, 6 e 7: iguais a Tabela 3 ate a 4a casa decimal;
%     Area 3: -409.33 MW (Tabela: -409.52); Area 8: -399.00 MW (-399.09);
%     Area 9 (referencia, nao controlada): 996.65 MW (994.59).
%   A barra 47 foi colocada na Area 3, que da o menor desvio.
%
%   Saida: area(k) = area da barra k (k = 1..118).


    %% a) Area de Cada Barra
    area = [ ...
        1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 ...
        2 2 4 4 4 4 4 4 4 2 2 4 3 3 3 3 3 3 3 3 ...
        3 3 9 9 9 9 3 9 9 9 9 8 8 8 8 8 8 8 8 8 ...
        8 8 8 9 9 9 9 9 9 4 4 4 4 5 5 5 5 5 7 7 ...
        9 5 5 5 5 5 5 5 7 7 7 7 7 7 7 7 7 7 7 6 ...
        6 6 6 6 6 6 6 6 6 6 6 6 2 4 4 9 1 5 ...
        ].';
end
