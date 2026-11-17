function S = Le_Series_ACTIVSg2000(pasta)
%LE_SERIES_ACTIVSG2000  Le as series horarias de 2016 do ACTIVSg2000 (Texas A&M).
%
%   S = Le_Series_ACTIVSg2000(pasta)
%
%   Arquivos esperados na pasta (nomes como baixados do site):
%     ACTIVISg2000_load_time_series_MW.csv
%     ACTIVISg2000_load_time_series_MVAR.csv
%     ACTIVISg2000_renewable_time_series_MW.csv
%
%   Saida (8784 horas):
%     S.data, S.hora     : data e hora de cada registro (texto)
%     S.barraCarga       : numero ORIGINAL das barras com carga (nC x 1)
%     S.Pd, S.Qd         : carga por barra [MW, Mvar] (8784 x nC, single)
%     S.barraRen         : numero original da barra de cada usina renovavel
%     S.tipoRen          : 'Wind' ou 'Solar' de cada usina
%     S.Pren             : geracao renovavel por usina [MW] (8784 x nR, single)
%     S.PdTotal          : carga total do sistema [MW]
%
%   Cargas da mesma barra (#1, #2, ...) sao somadas.
%
%   Referencia dos dados: H. Li, J. H. Yeo, A. L. Bornsheuer, T. J. Overbye,
%   "The Creation and Validation of Load Time Series for Synthetic Electric
%   Power Systems", IEEE Trans. Power Syst., vol. 36, no. 2, 2021.



    %% ____________________________________________________________________
    %%
    %%                           1. CARGAS
    %% ____________________________________________________________________

    %% a) Potencia Ativa
    [cab, dados, data, hora] = le_csv(fullfile(pasta, 'ACTIVISg2000_load_time_series_MW.csv'));

    % Colunas "Bus xxxx #y MW"
    colB  = find(strncmp(cab, 'Bus ', 4));
    bCol  = barra_do_cabecalho(cab(colB));
    [S.barraCarga, ~, j] = unique(bCol);
    nC    = numel(S.barraCarga);
    nh    = size(dados, 1);

    S.Pd = zeros(nh, nC, 'single');

    for i = 1:numel(colB)
        S.Pd(:, j(i)) = S.Pd(:, j(i)) + single(dados(:, colB(i)));
    end

    S.PdTotal = dados(:, 2);
    S.data    = data;
    S.hora    = hora;


    %% b) Potencia Reativa
    [cab, dados] = le_csv(fullfile(pasta, 'ACTIVISg2000_load_time_series_MVAR.csv'));
    colB = find(strncmp(cab, 'Bus ', 4));
    bCol = barra_do_cabecalho(cab(colB));
    [~, j] = ismember(bCol, S.barraCarga);

    S.Qd = zeros(nh, nC, 'single');

    for i = 1:numel(colB)
        S.Qd(:, j(i)) = S.Qd(:, j(i)) + single(dados(:, colB(i)));
    end



    %% ____________________________________________________________________
    %%
    %%                         2. RENOVAVEIS
    %% ____________________________________________________________________

    %% a) Leitura (a primeira linha traz o tipo de cada usina)
    arq = fullfile(pasta, 'ACTIVISg2000_renewable_time_series_MW.csv');
    fid = fopen(arq, 'r');
    l1  = fgetl(fid);
    fclose(fid);
    tipos = strsplit(l1, ',', 'CollapseDelimiters', false);

    [cab, dados] = le_csv(arq);

    % Colunas "Gen xxxx #y MW"
    colG = find(strncmp(cab, 'Gen ', 4));


    %% b) Estrutura de Saida
    S.barraRen = barra_do_cabecalho(cab(colG));
    S.tipoRen  = tipos(colG + 2);
    S.Pren     = single(dados(:, colG));
end



%% ________________________________________________________________________
%%
%%                       FUNCOES AUXILIARES
%% ________________________________________________________________________

function [cab, dados, data, hora] = le_csv(arq)
%LE_CSV  Le um CSV das series: 1a linha ignorada, 2a linha = cabecalho,
%        colunas 1-2 = data e hora (texto), demais = numeros.
%        dados(:, k) corresponde a coluna k+2 do arquivo.

    fid = fopen(arq, 'r');

    % Linha "PWOPFTimePoint,..."
    fgetl(fid);
    l2  = fgetl(fid);
    cab = strsplit(l2, ',', 'CollapseDelimiters', false);

    % Sem Date e Time
    cab = cab(3:end);
    nNum = numel(cab);

    fmt = ['%s %s', repmat(' %f', 1, nNum)];
    C   = textscan(fid, fmt, 'Delimiter', ',', 'EmptyValue', 0, 'CollectOutput', true);
    fclose(fid);

    data  = C{1}(:,1);
    hora  = C{1}(:,2);
    dados = C{2};
end


function b = barra_do_cabecalho(nomes)
%BARRA_DO_CABECALHO  Extrai o numero da barra de "Bus 1001 #1 MW" ou "Gen 1011 #1 MW".

    b = zeros(numel(nomes), 1);

    for i = 1:numel(nomes)
        p = strsplit(strtrim(nomes{i}), ' ');
        b(i) = str2double(p{2});
    end
end
