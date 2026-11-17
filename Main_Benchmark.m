%                                                                       JMJ
%
%       "Aquilo que Deus quer para nos e muito mais belo do que tudo
%              que poderiamos pedir com a nossa imaginacao."
%
%                                 Chiara Corbella Petrillo


%% ________________________________________________________________________
%%
%%                        1. INFORMACOES INICIAIS
%% ________________________________________________________________________

%% a) Comentarios Iniciais
% Este programa mede o tempo de processamento do laco de controle completo
% (todos os estagios, incluindo os reinicios) nos casos do artigo:
%   4 barras    : Casos 1, 2 e 3 de Main_4Barras;
%   IEEE 118    : caso base, transacao de +60 MW e exportacao de 250 MW
%                 com limite termico (troca de equacao), de Main_IEEE118;
%   ACTIVSg2000 : perda dos parques eolicos nas horas 5368 e 8504, de
%                 Main_ACTIVSg2000.
% Variantes:
%   HELM (Wynn)                      - reinicia da solucao germinal
%   NR esparso com partida plana     - reinicia do perfil plano
%   NR esparso com partida quente    - reinicia da solucao do estagio anterior
%
% Protocolo:
%   - Os dados (Ybus, areas, GCIs) sao montados fora da medicao;
%   - Cada caso e executado nrep vezes; reporta-se a mediana;
%   - Uma execucao de aquecimento e descartada (compilacao JIT do MATLAB);
%   - Mesma tolerancia (cfg.xi = cfg.tolNR) e mesmo laco de controle;
%   - Reporta-se tambem o custo: ordens (HELM) e iteracoes (NR), os
%     estagios e o tempo medio por estagio.
%
% Os tempos dependem da versao do MATLAB e do computador; o artigo reporta
% os resultados obtidos no MATLAB R2018a e no MATLAB R2026b.


%% b) Limpeza de Figuras e Variaveis
clc
clear all
close all
format short
warning off


%% c) Pastas das Funcoes e dos Dados
addpath('funcoes', 'dados');



%% ________________________________________________________________________
%%
%%                       2. CONFIGURACOES GLOBAIS
%% ________________________________________________________________________

%% a) Repeticoes
% Repeticoes por caso
nrep = 20;


%% b) Variantes Comparadas
var.nome   = {'HELM (Wynn)', 'NR plana', 'NR quente'};
var.metodo = {'HELM', 'NR', 'NR'};
var.part   = {'plana', 'plana', 'quente'};
nv = numel(var.nome);



%% ________________________________________________________________________
%%
%%                    3. PREPARACAO DOS CASOS (fora da medicao)
%% ________________________________________________________________________

%% a) Tutorial de 4 Barras
cfg = Config_CIA('liu4');
mpc = Carrega_Caso(cfg);

% [exportacao, limQ, limP]
par = [ 40, 0, 0;
        40, 1, 1;
       100, 1, 1];

nome  = {};
Dlist = {};

for c = 1:3
    cfgc = cfg;
    cfgc.NIesp_MW = -par(c,1);
    cfgc.limQ     = par(c,2) == 1;
    cfgc.limP     = par(c,3) == 1;
    Dlist{end+1}  = CIA_Dados(mpc, cfgc);
    nome{end+1}   = sprintf('4 barras, Caso %d', c);
end


%% b) IEEE 118 Barras com 3 Areas
cfg = Config_CIA('ieee118');
mpc = Carrega_Caso(cfg);
D0  = CIA_Dados(mpc, cfg);
tipo0 = D0.tipo;
tipo0(D0.gci) = 2;
D0.cfg.limQ = false;
V0  = CIA_Resolver(D0, 'NR', tipo0, D0.P0, D0.Q0, []);
NI0 = D0.baseMVA*CIA_Intercambio(D0, V0);

casos118.nome = {'IEEE 118, base', 'IEEE 118, +60 MW', 'IEEE 118, 250 MW (lim. termico)'};
casos118.NI   = [NI0(1),       NI0(2);
                 NI0(1) - 60,  NI0(2);
                 -250,         NI0(2)];
casos118.term = [false, false, true];

for c = 1:3
    cfgc = cfg;
    cfgc.NIesp_MW   = casos118.NI(c,:);
    cfgc.limTermico = casos118.term(c);
    Dlist{end+1}    = CIA_Dados(mpc, cfgc);
    nome{end+1}     = casos118.nome{c};
end


%% c) ACTIVSg2000: Perda dos Parques Eolicos
horas  = [5368 8504];
parque = [2127 4153];

for ih = 1:numel(horas)
    cfg      = Config_CIA('activsg_h');
    cfg.caso = sprintf('case_ACTIVSg2000_h%04d', horas(ih));
    mpc      = Carrega_Caso(cfg);
    wf       = zeros(1, numel(parque));

    for i = 1:numel(parque)
        wf(i) = find(mpc.numOrig == parque(i));
    end

    cfg0 = cfg;
    cfg0.NIesp_MW = zeros(1, 5);
    D0    = CIA_Dados(mpc, cfg0);
    tipo0 = D0.tipo;
    tipo0(D0.gci) = 2;
    V0    = CIA_Resolver(D0, 'NR', tipo0, D0.P0, D0.Q0, []);

    cfg.NIesp_MW = (D0.baseMVA*CIA_Intercambio(D0, V0)).';
    cfg.WF       = wf;
    cfg.WF_MW    = (D0.baseMVA*(D0.P0(wf) + D0.Pd(wf))).';
    Dlist{end+1} = CIA_Dados(mpc, cfg, [1, -1, -1]);
    nome{end+1}  = sprintf('ACTIVSg2000, hora %d', horas(ih));
end

nc = numel(Dlist);



%% ________________________________________________________________________
%%
%%                             4. MEDICAO
%% ________________________________________________________________________

%% a) Estrutura de Resultados
% Mediana do tempo [s]
tMed  = nan(nc, nv);

% Ordens ou iteracoes
custo = nan(nc, nv);

% Estagios
est   = nan(nc, nv);

% Solucao encontrada
ok    = false(nc, nv);

% Max|V_HELM - V_NR plana|
dV    = nan(nc, 1);


%% b) Laco de Medicao
for i = 1:nc
    Vref = [];

    for v = 1:nv
        D = Dlist{i};
        D.cfg.partidaNR = var.part{v};

        % Aquecimento (descartado)
        r  = CIA_Controle(D, var.metodo{v});
        tt = zeros(nrep, 1);

        for rep = 1:nrep
            t0 = tic;
            r  = CIA_Controle(D, var.metodo{v});
            tt(rep) = toc(t0);
        end

        tMed(i,v)  = median(tt);
        custo(i,v) = r.custo;
        est(i,v)   = r.estagios;
        ok(i,v)    = ~r.semSolucao;

        if v == 1
            Vref = r.V;
        elseif v == 2
            dV(i) = max(abs(Vref - r.V));
        end
    end

    fprintf('%-34s: %s s\n', nome{i}, mat2str(tMed(i,:), 3));
end



%% ________________________________________________________________________
%%
%%                             5. RESUMO
%% ________________________________________________________________________

%% a) Tabela
fprintf('\n%-34s %5s | %10s %10s %10s | %6s %6s %6s | %8s | %7s\n', 'Caso', 'est.', ...
    'HELM [ms]', 'NRp [ms]', 'NRq [ms]', 'ordens', 'itNRp', 'itNRq', 'NRp/HELM', 'max|dV|');

for i = 1:nc
    fprintf('%-34s %5d | %10.2f %10.2f %10.2f | %6d %6d %6d | %8.2f | %7.1e\n', nome{i}, est(i,1), ...
        1e3*tMed(i,1), 1e3*tMed(i,2), 1e3*tMed(i,3), custo(i,1), custo(i,2), custo(i,3), ...
        tMed(i,2)/tMed(i,1), dV(i));
end

if ~all(ok(:))
    fprintf('\nATENCAO: casos sem solucao em algum metodo: %s\n', mat2str(find(~all(ok,2)).'));
end


%% b) Tempo Medio por Estagio [ms]
fprintf('\n%-34s %12s %12s %12s\n', 'Caso', 'HELM/est.', 'NRp/est.', 'NRq/est.');

for i = 1:nc
    fprintf('%-34s %12.2f %12.2f %12.2f\n', nome{i}, 1e3*tMed(i,:)./est(i,:));
end


%% c) Gravacao
save('Resultados_Benchmark.mat', 'nome', 'tMed', 'custo', 'est', 'ok', 'dV', 'var', 'nrep');
