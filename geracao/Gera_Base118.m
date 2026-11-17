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
% Este programa gera o ponto de operacao base do IEEE 118 barras usado no
% estudo do CIA (arquivo case118_base.m). Requer o MATPOWER (runopf).
% O arquivo gerado ja acompanha o repositorio (pasta dados); este programa
% so precisa ser executado para mudar a reserva ou a margem termica.
%
% Procedimento:
%   1) Dados: pglib_opf_case118_ieee (IEEE PES PGLib-OPF v23.07), com a
%      mesma rede do case118 e capacidades de linha (rateA). Nele, os 35
%      compensadores sincronos tem Pmax = 0 (nao geram potencia ativa).
%   2) OPF AC com reserva para o CAG: os limites dos geradores reais sao
%      estreitados para [Pmin + r*Pmax, Pmax - r*Pmax], de modo que todo
%      gerador fique com pelo menos r*Pmax de margem para subir e para descer.
%   3) Margem termica: no OPF, a capacidade das linhas e reduzida para
%      (1 - rL)*rateA, de modo que o ponto base nao fique encostado no
%      limite termico (margem operativa para as transacoes).
%   4) O ponto de operacao (Pg, Vg, V, theta) e gravado com os limites
%      originais, que sao os usados pelo laco de saturacao do CIA e na
%      verificacao do carregamento.


%% b) Limpeza de Figuras e Variaveis
clc
clear all
close all
format short


%% c) Pasta dos Dados
% Os casos sao lidos e gravados na pasta dados do repositorio.
pastaDados = fullfile('..', 'dados');
addpath(pastaDados);



%% ________________________________________________________________________
%%
%%                       2. OPF COM RESERVA
%% ________________________________________________________________________

%% a) Dados e Reserva
% Reserva de 10% do Pmax
r   = 0.10;

% Margem termica de 10% da capacidade
rL  = 0.10;
mpc = pglib_opf_case118_ieee;

% Geradores reais (Pmax > 0)
ger = find(mpc.gen(:,9) > 0);


%% b) Limites Estreitados para o OPF
opf = mpc;
opf.gen(ger,9)  = mpc.gen(ger,9)  - r*mpc.gen(ger,9);
opf.gen(ger,10) = mpc.gen(ger,10) + r*mpc.gen(ger,9);
opf.branch(:,6:8) = (1 - rL)*mpc.branch(:,6:8);


%% c) Solucao do OPF
opcoes = mpoption('verbose', 0, 'out.all', 0);
res    = runopf(opf, opcoes);

if ~res.success
    error('Gera_Base118: o OPF nao convergiu.');
end



%% ________________________________________________________________________
%%
%%                     3. CASO BASE (LIMITES ORIGINAIS)
%% ________________________________________________________________________

%% a) Ponto de Operacao
base = mpc;
base.gen(:,2) = res.gen(:,2);
base.gen(:,3) = res.gen(:,3);
base.gen(:,6) = res.bus(res.gen(:,1), 8);
base.bus(:,8) = res.bus(:,8);
base.bus(:,9) = res.bus(:,9);


%% b) Verificacao
S = max(abs(res.branch(:,14) + 1j*res.branch(:,15)), abs(res.branch(:,16) + 1j*res.branch(:,17)));
fprintf('Custo: %.2f $/h | V em [%.3f, %.3f] p.u. | carregamento maximo: %.1f%%\n', ...
    res.f, min(res.bus(:,8)), max(res.bus(:,8)), 100*max(S./res.branch(:,6)));
fprintf('Geracao total: %.1f MW | carga: %.1f MW\n', sum(res.gen(:,2)), sum(res.bus(:,3)));


%% c) Gravacao do Arquivo de Dados
fid = fopen(fullfile(pastaDados, 'case118_base.m'), 'w');
fprintf(fid, 'function mpc = case118_base\n');
fprintf(fid, '%%CASE118_BASE  IEEE 118 barras (PGLib-OPF v23.07) no ponto base do estudo do CIA.\n');
fprintf(fid, '%%   Ponto de operacao obtido por OPF AC com reserva de %.0f%% do Pmax em\n', 100*r);
fprintf(fid, '%%   todos os geradores reais e margem termica de %.0f%% nas linhas\n', 100*rL);
fprintf(fid, '%%   (gerado por Gera_Base118.m). Limites originais.\n');
fprintf(fid, '%%   Fonte dos dados: IEEE PES Power Grid Library, licenca CC BY 4.0.\n\n');
fprintf(fid, 'mpc.version = ''2'';\n');
fprintf(fid, 'mpc.baseMVA = %g;\n\n', base.baseMVA);
Escreve_Matriz(fid, 'mpc.bus', base.bus);
Escreve_Matriz(fid, 'mpc.gen', base.gen);
Escreve_Matriz(fid, 'mpc.branch', base.branch);
fprintf(fid, 'end\n');
fclose(fid);

fprintf('Arquivo case118_base.m gravado.\n');

