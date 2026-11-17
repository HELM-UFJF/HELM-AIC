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
% Este programa reproduz os casos de Luo e Liu (EPSR 225, 109721, 2023) e
% os compara com a formulacao proposta:
%   Parte A: sistema de 4 barras (Tabelas 1 e 2 do artigo);
%   Parte B: IEEE 118 barras com 9 areas e intercambio programado igual a
%            0,8 x o original (Tabelas 3 e 4), sem e com limites de geracao;
%   Parte C: teste de estresse (Figs. 6 e 7): intercambio programado igual
%            a k x o original, sem e com limites de geracao.
%
% Convencao de sinal: no artigo a exportacao e positiva; no codigo, NI e a
% potencia recebida pela area. As tabelas impressas usam a convencao do
% artigo.


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
disp('.     AREA INTERCHANGE CONTROL: COMPARACAO COM LUO E LIU (2023)    .')
disp('.                                                                  .')
disp('....................................................................')



%% ________________________________________________________________________
%%
%%                   2. PARTE A: TUTORIAL DE 4 BARRAS
%% ________________________________________________________________________

%% a) Dados (Area 2 = barras 2 e 4; intercambio desejado: 40 MW)
cfg = Config_CIA('liu4');
mpc = Carrega_Caso(cfg);
D   = CIA_Dados(mpc, cfg);


%% b) Intercambio Original (sem CIA)
tipo0 = D.tipo;
tipo0(D.gci) = 2;
V0 = CIA_Resolver(D, 'NR', tipo0, D.P0, D.Q0, []);

fprintf('\n===== Parte A: 4 barras =====\n');
fprintf('Intercambio original da Area 2: %.4f MW  (Luo e Liu: 66.23 MW)\n', -D.baseMVA*CIA_Intercambio(D, V0));


%% c) Com CIA (HELM e NR)
metodos = {'HELM', 'NR'};

for m = 1:2
    r = CIA_Controle(D, metodos{m});
    fprintf('%-4s: intercambio %.4f MW  lambda %.4f MW  |V2| %.4f  |V3| %.4f  (%d ordens/iteracoes)\n', ...
        metodos{m}, -D.baseMVA*r.NI, D.baseMVA*r.lam, abs(r.V(2)), abs(r.V(3)), r.info.n);
end

fprintf('Luo e Liu: lambda -26.518 MW  |V2| 0.9825  |V3| 0.9691\n');



%% ________________________________________________________________________
%%
%%               3. PARTE B: IEEE 118 BARRAS COM 9 AREAS
%% ________________________________________________________________________

%% a) Dados
cfg    = Config_CIA('liu118');
mpc    = Carrega_Caso(cfg);

% Tabela 3 de Luo e Liu (2023) (exportacao)
NIorig = -cfg.NIesp_MW/0.8;
D      = CIA_Dados(mpc, cfg);


%% b) Intercambios Originais (sem CIA)
tipo0 = D.tipo;
tipo0(D.gci) = 2;
V0  = CIA_Resolver(D, 'NR', tipo0, D.P0, D.Q0, []);
NI0 = -D.baseMVA*CIA_Intercambio(D, V0);

fprintf('\n===== Parte B: IEEE 118, 9 areas =====\n');
fprintf('%-5s %14s %14s\n', 'Area', 'Calculado', 'Luo e Liu');

for ia = 1:D.nA
    fprintf('%-5d %14.4f %14.4f\n', D.areasCtrl(ia), NI0(ia), NIorig(ia));
end


%% c) Sem e Com Limites de Geracao
caso = {'sem limites', 'com limites'};

for c = 1:2
    D.cfg.limP = (c == 2);
    rH = CIA_Controle(D, 'HELM');
    rN = CIA_Controle(D, 'NR');
    RB.(sprintf('c%d', c)) = rH;

    fprintf('\n--- %s ---\n', caso{c});
    fprintf('HELM x NR: max|dV| = %.2e p.u.;  GCIs fixados: %d\n', max(abs(rH.V - rN.V)), numel(rH.fixados));

    % i. Fluxos nas linhas da Area 5 (Tabela 4 de Luo e Liu (2023))
    lig = D.lig(5);
    V   = rH.V;
    fprintf('  Linhas da Area 5 (exportacao, MW):\n');

    for i = 1:numel(lig.r)
        p = real(V(lig.r(i))*conj(lig.y(i)*(V(lig.r(i)) - V(lig.t(i))) + lig.ysh(i)*V(lig.r(i))));
        fprintf('    %3d-%3d: %10.4f\n', lig.r(i), lig.t(i), D.baseMVA*p);
    end

    % ii. Geradores fora dos limites (so faz sentido sem limites)
    if c == 1
        fprintf('  Barras de regulacao fora dos limites de geracao:\n');

        for k = D.gci.'
            if rH.Pg(k) < D.Pmin(k) - 1e-6 || rH.Pg(k) > D.Pmax(k) + 1e-6
                fprintf('    barra %3d (Area %d): Pg = %8.2f MW, limites [%g, %g] MW\n', k, D.area(k), ...
                    D.baseMVA*rH.Pg(k), D.baseMVA*D.Pmin(k), D.baseMVA*D.Pmax(k));
            end
        end
    end
end



%% ________________________________________________________________________
%%
%%                 4. PARTE C: TESTE DE ESTRESSE (k x original)
%% ________________________________________________________________________

%% a) Fatores
k  = 0.05:0.05:3.50;
nk = numel(k);

E.k   = k(:);

% Convergencia HELM (sem, com limites)
E.cH  = false(nk, 2);

% Convergencia NR
E.cN  = false(nk, 2);

% Maior residuo HELM [p.u.]
E.rH  = nan(nk, 2);

% Maior residuo NR [p.u.]
E.rN  = nan(nk, 2);

% Lambda da Area 5 (HELM) [MW]
E.lH  = nan(nk, 2);

% Lambda da Area 5 (NR) [MW]
E.lN  = nan(nk, 2);

% Areas sem controle (com limites)
E.loc = zeros(nk, 2);


%% b) Laco
for c = 1:2
    for i = 1:nk
        cfgk = cfg;
        cfgk.NIesp_MW = -k(i)*NIorig;
        cfgk.limP     = (c == 2);

        % Como em Luo e Liu (2023)
        cfgk.maxitNR  = 60;
        cfgk.nmax     = 60;
        Dk = CIA_Dados(mpc, cfgk);

        rH = CIA_Controle(Dk, 'HELM');
        rN = CIA_Controle(Dk, 'NR');

        E.cH(i,c)  = ~rH.semSolucao;
        E.cN(i,c)  = ~rN.semSolucao;
        E.rH(i,c)  = CIA_Residuo(Dk, rH);
        E.rN(i,c)  = CIA_Residuo(Dk, rN);
        E.lH(i,c)  = Dk.baseMVA*rH.lam(5);
        E.lN(i,c)  = Dk.baseMVA*rN.lam(5);
        E.loc(i,c) = sum(rN.loc);
    end

    fprintf('\nEstresse %s: ultimo k convergido HELM = %.2f, NR = %.2f\n', caso{c}, ...
        max([E.k(E.cH(:,c)); NaN]), max([E.k(E.cN(:,c)); NaN]));
end


%% c) Dados para o pgfplots
dat = [E.k, log10(E.rH), log10(E.rN), E.lH, E.lN, E.loc(:,2)];
fid = fopen('dados_Liu118_Estresse.dat', 'w');
fprintf(fid, 'k logResHELM_sem logResHELM_com logResNR_sem logResNR_com lamH_sem lamH_com lamN_sem lamN_com nLOC_com\n');
fprintf(fid, [repmat('%.6g ', 1, size(dat,2)-1), '%.6g\n'], dat.');
fclose(fid);

save('Resultados_Liu.mat', 'E', 'RB');
