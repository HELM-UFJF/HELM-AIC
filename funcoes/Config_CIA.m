function cfg = Config_CIA(versao)
%CONFIG_CIA  Configuracao dos estudos de Controle de Intercambio de Area (CIA).
%
%   cfg = Config_CIA(versao)
%
%   versao: 'liu4'      -> sistema de 4 barras de Luo e Liu (EPSR, 2023)
%           'liu118'    -> IEEE 118 barras com as 9 areas de Luo e Liu
%           'ieee118'   -> IEEE 118 barras (PGLib-OPF) com 3 areas
%           'activsg_h' -> ACTIVSg2000 em uma hora da serie de 2016
%
%   Todos os programas principais leem os parametros desta funcao. Os
%   intercambios especificados e as horas do ACTIVSg2000 sao definidos nos
%   proprios programas principais.



    %% ____________________________________________________________________
    %%
    %%                     1. PARAMETROS COMUNS
    %% ____________________________________________________________________

    %% a) Versao
    if nargin < 1 || isempty(versao)
        versao = 'ieee118';
    end

    cfg.versao = versao;


    %% b) Tolerancias e Limites de Ordem e de Iteracao
    % Tolerancia do HELM [p.u.]
    cfg.xi = 1e-6;

    % Tolerancia de desbalanco do NR [p.u.]
    cfg.tolNR = 1e-6;

    % Ordem maxima da serie
    cfg.nmax = 60;

    % Criterio de parada do HELM:
    %   'passo'  : diferenca entre estimativas sucessivas
    %   'residuo': desbalanco de P, Q e |V| (mesmo criterio do NR)
    cfg.critHELM = 'passo';

    % Avaliacao das series: 'Wynn' (incremental) ou 'WynnCompleto'
    cfg.serie = 'Wynn';

    % Parada por estagnacao (ordens sem melhora)
    cfg.nestag = 15;

    % Iteracoes maximas do NR por estagio
    cfg.maxitNR = 30;

    % Partida do NR em cada estagio: 'plana' ou 'quente'
    cfg.partidaNR = 'plana';


    %% c) Laco Externo de Controle
    % Violacoes minimas consideradas [p.u.]
    cfg.epsP = 1e-9;
    cfg.epsQ = 1e-4;

    % Numero maximo de estagios
    cfg.maxEst = 300;

    % Ordem de verificacao dos limites: 'QP' (Q antes de P) ou 'PQ'
    cfg.ordemLim = 'QP';

    % Barras com Q fixado por estagio (Inf = todas as violacoes)
    cfg.nQporEstagio = Inf;

    % Retorno PQ -> PV: folga de tensao [p.u.] e numero maximo por barra
    cfg.retornoQ    = true;
    cfg.tolVret     = 1e-4;
    cfg.maxRetornoQ = 2;

    % Limites de Q e de P dos geradores
    cfg.limQ = true;
    cfg.limP = true;

    % Limites de P: 'um' (um GCI por estagio) ou 'todos' (como em Luo e Liu)
    cfg.estrategia = 'um';

    % Limite termico no CIA (troca de equacao), tolerancia e numero
    % maximo de ajustes ou trocas por area
    cfg.limTermico = false;
    cfg.tolTerm    = 1e-4;
    cfg.maxTerm    = 10;


    %% d) Areas e Fatores de Participacao
    % Papel da area da barra de referencia (ver CIA_Beta):
    %   'folga'  : a area regula a propria geracao (CAG proprio)
    %   'nenhum' : os GCIs da area nao participam
    %   'espelho': os GCIs da area compensam as demais areas
    cfg.modoArea1 = 'folga';

    % Fatores de participacao iguais na area (true) ou proporcionais ao
    % despacho (false)
    cfg.betaIgual = false;

    % Geracao programada da barra de referencia [MW] ([] = a do caso)
    cfg.PrefEsp_MW = [];

    % Intercambio liquido recebido por cada area controlada [MW], na ordem
    % crescente do numero da area (negativo = a area exporta)
    cfg.NIesp_MW = [];


    %% e) Caso, Areas e GCIs
    cfg.caso      = '';
    cfg.mapaAreas = [];
    cfg.combGCI   = {};
    cfg.GCIextra  = [];


    %% f) Parques Eolicos (indices das barras e potencia nominal [MW])
    cfg.WF    = [];
    cfg.WF_MW = [];



    %% ____________________________________________________________________
    %%
    %%                     2. PARAMETROS POR VERSAO
    %% ____________________________________________________________________

    switch lower(versao)

        %% a) Sistema de 4 Barras de Luo e Liu (EPSR, 2023)
        % Area 1 (referencia) = barras 1 e 3; Area 2 = barras 2 e 4.
        % GCI: gerador da barra 4. Exportacao programada da Area 2: 40 MW.
        case 'liu4'
            cfg.caso       = 'case4gs';
            cfg.estrategia = 'todos';
            cfg.limQ       = false;
            cfg.mapaAreas  = [1 2 1 2];
            cfg.GCIextra   = 4;
            cfg.NIesp_MW   = -40;
            cfg.modoArea1  = 'nenhum';

        %% b) IEEE 118 Barras com as 9 Areas de Luo e Liu (EPSR, 2023)
        % Areas reconstruidas da Fig. 4 do artigo (Areas118_Liu), com a
        % Area 9 como referencia. GCIs: barras de regulacao do artigo, com
        % fatores iguais. Intercambio programado = 0,8 x o original.
        case 'liu118'
            cfg.caso       = 'case118';
            cfg.estrategia = 'todos';
            cfg.limQ       = false;
            cfg.mapaAreas  = Areas118_Liu();
            cfg.GCIextra   = [1 4 6 8 10 12, 15 18 19 31 113, 34 36 40 42, ...
                              24 25 26 27 32 70 72 73, 74 76 77 85 87, ...
                              100 103 104 105 107 110 111 112, 80 89 90 91 92 99, ...
                              54 55 56 59 61 62];

            % Exportacao original de cada area (Tabela 3 de Luo e Liu) [MW]
            NIorig = [135.0676 -366.5963 -409.5160 215.0399 -528.7555 ...
                      -50.9500 448.6331 -399.0895];

            cfg.NIesp_MW  = -0.8*NIorig;
            cfg.modoArea1 = 'nenhum';
            cfg.betaIgual = true;

        %% c) IEEE 118 Barras (PGLib-OPF) com 3 Areas
        % Zonas do arquivo original (Areas118_3); a Zona 3 contem a barra
        % de referencia 69. GCIs: somente unidades com Pmax > 0.
        % Ponto base: OPF AC com reserva (case118_base, Gera_Base118).
        case 'ieee118'
            cfg.caso      = 'case118_base';
            cfg.mapaAreas = Areas118_3();
            cfg.GCIextra  = [10 12 25 26 31, 46 49 54 59 61 65 66, 80 87 89 100 103 111];
            cfg.NIesp_MW  = [0 0];

        %% d) ACTIVSg2000 em uma Hora da Serie de 2016
        % Areas originais 1 a 3 agrupadas na Area 1; demais mantidas (6
        % areas; a barra de referencia 7098 fica na Area 5). GCIs: unidades
        % a gas e a carvao. Ponto base: OPF AC (Gera_BaseACTIVSg2000).
        case 'activsg_h'
            cfg.caso       = 'case_ACTIVSg2000_h5368';
            cfg.estrategia = 'todos';
            cfg.mapaAreas  = [1 1 1 2 3 4 5 6];
            cfg.combGCI    = {'ng', 'coal'};
            cfg.limTermico = true;

        otherwise
            error('Config_CIA: versao "%s" desconhecida.', versao);
    end
end
