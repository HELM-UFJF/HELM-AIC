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
% Este programa gera os pontos de operacao do ACTIVSg2000 usados no estudo
% do CIA (arquivos case_ACTIVSg2000_hNNNN.m). Requer o MATPOWER (runopf) e
% as series horarias de 2016 (Series_ACTIVSg2000.mat, Le_Series_ACTIVSg2000).
% Os arquivos gerados ja acompanham o repositorio (pasta dados).
%
% Procedimento para cada hora escolhida:
%   1) Carga: cada uma das 8 areas originais recebe o total da serie na
%      hora, f_a(h) = soma da serie na area / soma do caso na area, com a
%      distribuicao e o fator de potencia do caso (Pd e Qd escalados por
%      f_a). Nao se usa a serie barra a barra: nas areas 5 e 6 os valores
%      por barra nao tem correlacao com as cargas do caso (coeficiente
%      0,00 e 0,06 no pico, contra 0,97 a 0,99 nas demais), embora os
%      totais coincidam, o que indica numeracao diferente nessas areas;
%   2) Eolica: os parques das barras 2127 (Oeste) e 4153 (area 4) recebem a
%      previsao da hora, CF(h) x potencia nominal do caso; CF(h) = serie
%      normalizada pelo seu maximo anual. As demais renovaveis ficam no
%      despacho do caso (apenas 2 parques variam, por escolha do estudo);
%   3) Compensacao reativa: em carga leve (f_a < 0,90) os bancos de
%      capacitores da area sao desligados, como na operacao real. Sem isso
%      o OPF da hora 8504 e inviavel: o carregamento das linhas gera mais
%      reativo do que os geradores conseguem absorver (Qmin);
%   4) Despacho: OPF AC com reserva de 10% do Pmax nas unidades que
%      participam do CIA (gas e carvao; as hidro do caso tem Pmin = Pmax)
%      e margem termica de 5% nas linhas (com 10% o OPF do pico ja e
%      inviavel: o caso original opera perto das capacidades);
%   5) Gravacao com os limites originais (o CIA e verificado contra eles).
%
% Horas (serie de 2016):
%   5368 - 11/08 15h: pico anual de carga (66,3 GW), eolica moderada;
%   8504 - 20/12 07h: carga de 70% do pico com os dois parques no maximo.
% Horas mais leves (abaixo de ~60% do pico) nao foram usadas: com todas as
% unidades do caso ligadas, a soma dos minimos tecnicos fica proxima da
% carga, e essas horas exigiriam um compromisso de unidades, fora do
% escopo deste estudo.


%% b) Limpeza de Figuras e Variaveis
clc
clear all
close all
format short
warning off


%% c) Pasta dos Dados
% Os casos sao lidos e gravados na pasta dados do repositorio.
pastaDados = fullfile('..', 'dados');
addpath(pastaDados);



%% ________________________________________________________________________
%%
%%                       2. DADOS E PARAMETROS
%% ________________________________________________________________________

%% a) Parametros
horas  = [5368 8504];

% Reserva de 10% do Pmax
r      = 0.10;

% Margem termica de 5% da capacidade
rL     = 0.05;

% Abaixo disso, capacitores desligados
fCap   = 0.90;

% Barras dos parques eolicos
parque = [2127 4153];


%% b) Caso e Series
m = case_ACTIVSg2000;
S = load('Series_ACTIVSg2000.mat');
on   = m.gen(:,8) == 1;
eol  = strcmp(m.genfuel(:), 'wind');

% Unidades do CIA
conv = on & ismember(m.genfuel(:), {'ng', 'coal'}) & m.gen(:,9) > 0;



%% ________________________________________________________________________
%%
%%                       3. OPF PARA CADA HORA
%% ________________________________________________________________________

for h = horas

    %% a) Carga da Hora (total da serie por area original)
    opf = m;
    [~, lin] = ismember(S.barraCarga, m.bus(:,1));
    aSer = m.bus(lin, 7);
    fA   = zeros(1, 8);

    for a = 1:8
        naArea = m.bus(:,7) == a;
        fA(a)  = sum(double(S.Pd(h, aSer == a)))/sum(m.bus(naArea,3));
        opf.bus(naArea,3) = fA(a)*m.bus(naArea,3);
        opf.bus(naArea,4) = fA(a)*m.bus(naArea,4);

        if fA(a) < fCap
            % Capacitores desligados
            opf.bus(naArea & m.bus(:,6) > 0, 6) = 0;
        end
    end


    %% b) Previsao dos Parques Eolicos
    Pprev = zeros(1, numel(parque));

    for i = 1:numel(parque)
        u  = S.barraRen == parque(i);
        w  = sum(S.Pren(:,u), 2);
        CF = w(h)/max(w);
        g  = find(m.gen(:,1) == parque(i) & eol & on);
        Pprev(i) = CF*sum(m.gen(g,9));
        opf.gen(g,2)  = CF*m.gen(g,9);

        % Previsao fixa no OPF
        opf.gen(g,9)  = CF*m.gen(g,9);
        opf.gen(g,10) = CF*m.gen(g,9);
    end


    %% c) Reserva e Margem Termica
    opf.gen(conv,9)  = m.gen(conv,9) - r*m.gen(conv,9);
    opf.gen(conv,10) = min(m.gen(conv,10) + r*m.gen(conv,9), opf.gen(conv,9));
    opf.branch(:,6:8) = (1 - rL)*m.branch(:,6:8);


    %% d) Solucao do OPF
    opcoes = mpoption('verbose', 0, 'out.all', 0);
    res    = runopf(opf, opcoes);

    if ~res.success
        error('Gera_BaseACTIVSg2000: o OPF da hora %d nao convergiu.', h);
    end


    %% e) Ponto de Operacao com os Limites Originais
    base = m;
    base.bus(:,3)  = opf.bus(:,3);
    base.bus(:,4)  = opf.bus(:,4);
    base.bus(:,6)  = opf.bus(:,6);
    base.bus(:,8)  = res.bus(:,8);
    base.bus(:,9)  = res.bus(:,9);
    base.gen(:,2)  = res.gen(:,2);
    base.gen(:,3)  = res.gen(:,3);
    [~, ig] = ismember(m.gen(:,1), m.bus(:,1));

    % Tensao especificada = OPF
    base.gen(:,6)  = res.bus(ig, 8);


    %% f) Verificacao
    Sr = max(abs(res.branch(:,14) + 1j*res.branch(:,15)), abs(res.branch(:,16) + 1j*res.branch(:,17)));
    fprintf('Hora %d: f_a = %s\n', h, mat2str(fA, 3));
    fprintf('Hora %d: carga %.0f MW | eolica (2 parques) [%.1f %.1f] MW | custo %.0f $/h\n', ...
        h, sum(opf.bus(:,3)), Pprev, res.f);
    fprintf('          V em [%.3f, %.3f] p.u. | carregamento maximo %.1f%%\n', ...
        min(res.bus(:,8)), max(res.bus(:,8)), 100*max(Sr./m.branch(:,6)));


    %% g) Gravacao
    nome = sprintf('case_ACTIVSg2000_h%04d', h);
    fid  = fopen(fullfile(pastaDados, [nome, '.m']), 'w');
    fprintf(fid, 'function mpc = %s\n', nome);
    fprintf(fid, '%%%s  ACTIVSg2000 na hora %d da serie de 2016 (estudo do CIA).\n', upper(nome), h);
    fprintf(fid, '%%   Carga da serie horaria; parques 2127 e 4153 na previsao da hora;\n');
    fprintf(fid, '%%   despacho por OPF AC com reserva de %.0f%% e margem termica de %.0f%%\n', 100*r, 100*rL);
    fprintf(fid, '%%   (gerado por Gera_BaseACTIVSg2000.m). Limites originais.\n');
    fprintf(fid, '%%   Fonte: Texas A&M ACTIVSg2000 (MATPOWER), licenca CC BY 4.0.\n\n');
    fprintf(fid, 'mpc.version = ''2'';\n');
    fprintf(fid, 'mpc.baseMVA = %g;\n\n', base.baseMVA);
    Escreve_Matriz(fid, 'mpc.bus', base.bus);
    Escreve_Matriz(fid, 'mpc.gen', base.gen);
    Escreve_Matriz(fid, 'mpc.branch', base.branch);
    fprintf(fid, 'mpc.genfuel = {\n');
    fprintf(fid, '\t''%s'';\n', base.genfuel{:});
    fprintf(fid, '};\n');
    fprintf(fid, 'end\n');
    fclose(fid);

    fprintf('          Arquivo %s.m gravado.\n', nome);
end
