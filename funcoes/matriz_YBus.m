function [YBus, G, B] = matriz_YBus(baseMVA, bus, gen, branch, FE, FC)
%MATRIZ_YBUS  Monta a matriz de admitancia nodal em formato esparso.
%
%   A matriz e montada de forma vetorizada com sparse(), o que reduz o custo
%   de montagem e permite usar fatoracoes esparsas no HELM e no NR.
%
%   Entradas: baseMVA, bus, gen, branch (formato de Converte_MATPOWER), FE, FC
%   Saidas  : YBus (esparsa, p.u.), G = real(YBus), B = imag(YBus)


    %% a) Leitura dos Dados
    [De, Para, ~, ~, ~, Y, ~, SHL, a, SHB, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, Nbarra, ~] = ...
        ler_dados(baseMVA, bus, gen, branch, FE, FC);


    %% b) Elementos da Diagonal
    % Barra De : Y/a^2 + shunt da linha
    % Barra Para: Y     + shunt da linha
    diagDe   = Y./(a.^2) + SHL;
    diagPara = Y + SHL;


    %% c) Elementos Fora da Diagonal
    fora = -Y./a;


    %% d) Montagem Esparsa (entradas repetidas sao somadas por sparse)
    lin = [De; Para; De; Para; (1:Nbarra).'];
    col = [De; Para; Para; De; (1:Nbarra).'];
    val = [diagDe; diagPara; fora; fora; SHB];

    YBus = sparse(lin, col, val, Nbarra, Nbarra);


    %% e) Condutancia e Susceptancia
    G = real(YBus);
    B = imag(YBus);
end
