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
% Este programa executa o estudo do CIA no IEEE 118 barras com as 3 zonas
% do arquivo de dados original (Areas118_3):
%   - GCIs: somente geradores reais (compensadores sincronos nao regulam P);
%   - fatores de participacao proporcionais a geracao nominal;
%   - limites de P dos GCIs e de Q de todos os geradores (com retorno PQ->PV);
%   - a area da folga (Zona 3) regula a propria geracao, mantendo a barra 69
%     no despacho do caso (modo 'folga').
%
%   Parte A: caso base (intercambios originais) e transacao programada em
%            que a Area 1 vende 60 MW a mais;
%   Parte B: varredura da exportacao da Area 1 (transferencia crescente),
%            com verificacao fisica: tensoes [0.94, 1.06] p.u., carregamento
%            dos ramos (capacidades do PGLib-OPF v23.07) e limites dos GCIs;
%   Parte C: mesma varredura com o limite termico no CIA: quando um ramo
%            ultrapassa a capacidade, a equacao de intercambio da area e
%            trocada pela equacao de fluxo do ramo no limite, e a
%            exportacao alcancada passa a ser o limite de transferencia.
%
% Dados: case118_base = IEEE 118 do PGLib-OPF v23.07 (mesma rede do case118,
% com capacidades de linha e compensadores com Pmax = 0) no ponto base
% obtido por OPF AC com 10% de reserva nos geradores e 10% de margem
% termica nos ramos (Gera_Base118.m).


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
disp('.        AREA INTERCHANGE CONTROL: IEEE 118 BARRAS, 3 AREAS        .')
disp('.                                                                  .')
disp('....................................................................')



%% ________________________________________________________________________
%%
%%                       2. CONFIGURACOES GLOBAIS
%% ________________________________________________________________________

%% a) Sistema e Configuracao
cfg   = Config_CIA('ieee118');
mpc   = Carrega_Caso(cfg);

% Capacidade dos ramos [MVA] (PGLib)
rateA = mpc.rateA;


%% b) Intercambios Originais (sem CIA, sem limites de Q)
D0 = CIA_Dados(mpc, cfg);
tipo0 = D0.tipo;
tipo0(D0.gci) = 2;
D0.cfg.limQ = false;
V0  = CIA_Resolver(D0, 'NR', tipo0, D0.P0, D0.Q0, []);

% Recebido pelas Areas 1 e 2 [MW]
NI0 = D0.baseMVA*CIA_Intercambio(D0, V0);

fprintf('\nIntercambios originais (recebidos): Area 1 = %.2f MW, Area 2 = %.2f MW\n', NI0(1), NI0(2));



%% ________________________________________________________________________
%%
%%                    3. PARTE A: CASO BASE E TRANSACAO
%% ________________________________________________________________________

%% a) Casos
casos.nome = {'Base (intercambios originais)', 'Transacao: Area 1 vende +60 MW'};
casos.NI   = [NI0(1),       NI0(2);
              NI0(1) - 60,  NI0(2)];


%% b) Laco nos Casos
for c = 1:2
    cfgc = cfg;
    cfgc.NIesp_MW = casos.NI(c,:);
    D = CIA_Dados(mpc, cfgc);

    rH = CIA_Controle(D, 'HELM');
    rN = CIA_Controle(D, 'NR');
    cg = CIA_Carregamento(mpc, rH.V, rateA);
    Vm = abs(rH.V);

    fprintf('\n===== %s =====\n', casos.nome{c});
    fprintf('HELM: %d estagios, %d ordens | NR: %d iteracoes | max|dV| = %.2e p.u.\n', ...
        rH.estagios, rH.custo, rN.custo, max(abs(rH.V - rN.V)));
    fprintf('NI = [%.3f %.3f] MW (especificado [%.3f %.3f])\n', D.baseMVA*rH.NI, casos.NI(c,:));
    fprintf('lambda: Area 1 = %.3f MW, Area 2 = %.3f MW, Area 3 (folga) = %.3f MW\n', ...
        D.baseMVA*rH.lam, D.baseMVA*rH.lamRef);
    fprintf('Folga: %.2f MW | Q no limite: %d barras | V em [%.3f, %.3f] p.u. | carregamento max: %.1f%%\n', ...
        D.baseMVA*rH.Pslack, numel(rH.fixQ), min(Vm), max(Vm), 100*max(cg));

    fprintf('  %-5s %-5s %12s %12s %10s\n', 'GCI', 'Area', 'Pg nom [MW]', 'Pg fin [MW]', 'Pmax [MW]');

    for k = D.gci.'
        fprintf('  %-5d %-5d %12.2f %12.2f %10.1f\n', k, D.area(k), D.baseMVA*D.Pg_nom(k), ...
            D.baseMVA*rH.Pg(k), D.baseMVA*D.Pmax(k));
    end

    RA.(sprintf('c%d', c)) = rH;
end



%% ________________________________________________________________________
%%
%%               4. PARTE B: TRANSFERENCIA CRESCENTE DA AREA 1
%% ________________________________________________________________________

%% a) Exportacao Especificada da Area 1 [MW]
exp1 = 0:25:800;
ne   = numel(exp1);

T.exp1   = exp1(:);

% HELM convergiu
T.cH     = false(ne,1);

% NR convergiu
T.cN     = false(ne,1);

% Area 1 sem controle
T.loc1   = false(ne,1);

% GCIs fixados
T.nsat   = nan(ne,1);

% Exportacao alcancada [MW]
T.NI1    = nan(ne,1);
T.Vmin   = nan(ne,1);
T.Vmax   = nan(ne,1);

% Maior carregamento [p.u.]
T.carg   = nan(ne,1);

% Ramo mais carregado
T.ramo   = nan(ne,2);
T.dV     = nan(ne,1);


%% b) Laco
for i = 1:ne
    cfgi = cfg;
    cfgi.NIesp_MW = [-exp1(i), NI0(2)];
    D = CIA_Dados(mpc, cfgi);

    rH = CIA_Controle(D, 'HELM');
    rN = CIA_Controle(D, 'NR');

    T.cH(i)   = ~rH.semSolucao;
    T.cN(i)   = ~rN.semSolucao;

    if T.cN(i)
        cg = CIA_Carregamento(mpc, rN.V, rateA);
        [T.carg(i), L] = max(cg);
        T.ramo(i,:) = mpc.branch(L, 1:2);
        T.loc1(i)  = rN.loc(1);
        T.nsat(i)  = numel(rN.fixados);
        T.NI1(i)   = -D.baseMVA*rN.NI(1);
        T.Vmin(i)  = min(abs(rN.V));
        T.Vmax(i)  = max(abs(rN.V));
    end

    if T.cH(i) && T.cN(i)
        T.dV(i) = max(abs(rH.V - rN.V));
    end

    fprintf('exp %4d MW: HELM %d NR %d | alcancado %7.2f MW | LOC1 %d | fixados %2d | V [%.3f %.3f] | carg %5.1f%% (%d-%d) | dV %.1e\n', ...
        exp1(i), T.cH(i), T.cN(i), T.NI1(i), T.loc1(i), T.nsat(i), T.Vmin(i), T.Vmax(i), 100*T.carg(i), T.ramo(i,1), T.ramo(i,2), T.dV(i));
end


%% c) Limites Fisicos da Transferencia
iV = find(T.cN & (T.Vmin < 0.94 - 1e-3 | T.Vmax > 1.06 + 1e-3), 1);
iC = find(T.cN & T.carg > 1 + 1e-3, 1);
iL = find(T.cN & T.loc1, 1);

fprintf('\nPrimeira exportacao com tensao fora de [0.94, 1.06]: %s MW\n', num2str(exp1(iV)));
fprintf('Primeira exportacao com ramo acima de 100%%: %s MW\n', num2str(exp1(iC)));
fprintf('Primeira exportacao com perda de controlabilidade da Area 1: %s MW\n', num2str(exp1(iL)));


%% d) Dados para o pgfplots
dat = [T.exp1, T.NI1, double(T.loc1), T.nsat, T.Vmin, T.Vmax, 100*T.carg];
fid = fopen('dados_IEEE118_Transferencia.dat', 'w');
fprintf(fid, 'exp alcancado loc1 nsat Vmin Vmax carg\n');
fprintf(fid, [repmat('%.6g ', 1, size(dat,2)-1), '%.6g\n'], dat.');
fclose(fid);



%% ________________________________________________________________________
%%
%%          5. PARTE C: TRANSFERENCIA COM LIMITE TERMICO NO CIA
%% ________________________________________________________________________

%% a) Estrutura dos Resultados
TC.exp1  = exp1(:);
TC.cH    = false(ne,1);
TC.cN    = false(ne,1);

% Equacao do ramo ativa na Area 1
TC.ativo = false(ne,1);

% Ramo que limita a Area 1
TC.ramo  = nan(ne,2);

% Exportacao alcancada [MW]
TC.NI1   = nan(ne,1);

% Maior carregamento [p.u.]
TC.carg  = nan(ne,1);
TC.Vmin  = nan(ne,1);
TC.Vmax  = nan(ne,1);
TC.nsat  = nan(ne,1);
TC.dV    = nan(ne,1);

% |NI HELM - NI NR| [MW]
TC.dNI   = nan(ne,1);


%% b) Laco
for i = 1:ne
    cfgi = cfg;
    cfgi.NIesp_MW   = [-exp1(i), NI0(2)];
    cfgi.limTermico = true;
    D = CIA_Dados(mpc, cfgi);

    rH = CIA_Controle(D, 'HELM');
    rN = CIA_Controle(D, 'NR');

    TC.cH(i) = ~rH.semSolucao;
    TC.cN(i) = ~rN.semSolucao;

    if TC.cN(i)
        [~, ~, cg] = CIA_FluxoRamos(D, rN.V);
        TC.carg(i)  = max(cg);
        TC.ativo(i) = rN.termL(1) > 0;
        TC.NI1(i)   = -D.baseMVA*rN.NI(1);
        TC.Vmin(i)  = min(abs(rN.V));
        TC.Vmax(i)  = max(abs(rN.V));
        TC.nsat(i)  = numel(rN.fixados);

        if TC.ativo(i)
            L = rN.termL(1);
            TC.ramo(i,:) = [D.ramo.de(L), D.ramo.para(L)];
        end
    end

    if TC.cH(i) && TC.cN(i)
        TC.dV(i)  = max(abs(rH.V - rN.V));
        TC.dNI(i) = D.baseMVA*max(abs(rH.NI - rN.NI));
    end

    fprintf('exp %4d MW: HELM %d NR %d | alcancado %7.2f MW | ramo %s | fixados %2d | V [%.3f %.3f] | carg %6.2f%% | dV %.1e\n', ...
        exp1(i), TC.cH(i), TC.cN(i), TC.NI1(i), mat2str(TC.ramo(i,:)), TC.nsat(i), ...
        TC.Vmin(i), TC.Vmax(i), 100*TC.carg(i), TC.dV(i));

    for k = 1:numel(rN.infoT)
        fprintf('    aviso: %s\n', rN.infoT{k});
    end
end


%% c) Limite de Transferencia
iT = find(TC.ativo, 1);

if ~isempty(iT)
    fprintf('\nLimite de transferencia da Area 1 (ramo %d-%d no limite): %.2f MW\n', ...
        TC.ramo(iT,1), TC.ramo(iT,2), TC.NI1(iT));
end


%% d) Dados para o pgfplots
dat = [TC.exp1, TC.NI1, double(TC.ativo), TC.nsat, TC.Vmin, TC.Vmax, 100*TC.carg];
fid = fopen('dados_IEEE118_LimiteTermico.dat', 'w');
fprintf(fid, 'exp alcancado ativo nsat Vmin Vmax carg\n');
fprintf(fid, [repmat('%.6g ', 1, size(dat,2)-1), '%.6g\n'], dat.');
fclose(fid);

save('Resultados_IEEE118.mat', 'RA', 'T', 'TC', 'NI0');
