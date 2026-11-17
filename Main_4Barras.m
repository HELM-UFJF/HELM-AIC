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
% Este programa executa o CIA no sistema de 4 barras de Grainger e
% Stevenson, o mesmo da Secao 4.1 de Luo e Liu (EPSR, 2023). A Area 1
% (referencia) contem as barras 1 e 3; a Area 2 contem as barras 2 e 4,
% com o GCI da barra 4 (Pg = Pmax = 318 MW).
%   Caso 1: exportacao de 40 MW pela Area 2, sem limites (validacao);
%   Caso 2: o mesmo intercambio com os limites de Q dos geradores do caso
%           (|Qg| <= 100 Mvar); na solucao sem limites, Qg4 = 185 Mvar;
%   Caso 3: exportacao de 100 MW com os limites de P e Q. O GCI ja opera
%           em Pmax e nao pode aumentar a geracao: a Area 2 perde a
%           controlabilidade e o intercambio alcancado e informado.


%% b) Limpeza de Figuras e Variaveis
clc
clear all
close all
format short
warning off


%% c) Pastas das Funcoes e dos Dados
addpath('funcoes', 'dados');


%% d) Cabecalho
disp('....................................................................')
disp('.                                                                  .')
disp('.            AREA INTERCHANGE CONTROL: TUTORIAL 4 BARRAS           .')
disp('.                                                                  .')
disp('....................................................................')



%% ________________________________________________________________________
%%
%%                        2. DADOS E CASOS
%% ________________________________________________________________________

%% a) Sistema
cfg = Config_CIA('liu4');
mpc = Carrega_Caso(cfg);
metodos = {'HELM', 'NR'};


%% b) Casos: [exportacao (MW), limite de Q, limite de P]
casos.nome = {'Caso 1: 40 MW sem limites (Luo e Liu)', ...
              'Caso 2: 40 MW com limites de Q', ...
              'Caso 3: 100 MW com limites de P e Q'};
casos.par  = [ 40, 0, 0;
               40, 1, 1;
              100, 1, 1];


%% c) Intercambio Original (sem CIA)
D     = CIA_Dados(mpc, cfg);
tipo0 = D.tipo;
tipo0(D.gci) = 2;
V0 = CIA_Resolver(D, 'NR', tipo0, D.P0, D.Q0, []);

fprintf('\nExportacao original da Area 2: %.3f MW (Luo e Liu: 66.23 MW)\n', -D.baseMVA*CIA_Intercambio(D, V0));



%% ________________________________________________________________________
%%
%%                          3. SIMULACOES
%% ________________________________________________________________________

for c = 1:size(casos.par, 1)

    %% a) Configuracao do Caso
    cfgc = cfg;
    cfgc.NIesp_MW = -casos.par(c,1);
    cfgc.limQ     = casos.par(c,2) == 1;
    cfgc.limP     = casos.par(c,3) == 1;
    D = CIA_Dados(mpc, cfgc);

    fprintf('\n===== %s =====\n', casos.nome{c});


    %% b) HELM e NR
    for m = 1:2
        r  = CIA_Controle(D, metodos{m});
        Vm = abs(r.V);

        fprintf('%-4s: %d estagios | exportacao %.3f MW | lambda %.3f MW | LOC %d\n', metodos{m}, ...
            r.estagios, -D.baseMVA*r.NI, D.baseMVA*r.lam, r.loc);
        fprintf('      Pg4 %.2f MW | Qg4 %.2f Mvar | Qg1 %.2f Mvar | |V| = %s p.u.\n', ...
            D.baseMVA*r.Pg(4), D.baseMVA*r.Qg(4), D.baseMVA*r.Qg(1), mat2str(Vm.', 4));

        R(c).(metodos{m}) = r;
    end

    fprintf('max|V_HELM - V_NR| = %.2e p.u.\n', max(abs(R(c).HELM.V - R(c).NR.V)));
end

fprintf('\nLuo e Liu: lambda -26.518 MW | |V2| 0.9825 | |V3| 0.9691\n');
save('Resultados_4Barras.mat', 'R');
