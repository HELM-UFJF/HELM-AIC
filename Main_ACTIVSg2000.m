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
% Este programa executa o estudo do CIA no ACTIVSg2000 (2000 barras, Texas
% sintetico) em duas horas da serie de 2016:
%   5368 - 11/08 15h: pico anual de carga (66,3 GW);
%   8504 - 20/12 07h: 70% do pico, com os dois parques eolicos no maximo.
% Ponto base de cada hora: OPF AC com reserva (Gera_BaseACTIVSg2000.m).
% Intercambios programados = intercambios do ponto base.
%
% Evento: perda total da geracao dos parques eolicos das barras 2127
% (Area 1) e 4153 (Area 2), como no corte por vento forte. Compara-se:
%   - sem CIA: a folga (Area 5) assume a perda e os intercambios desviam;
%   - com CIA (HELM e NR): cada area repoe a propria perda com seus GCIs
%     (gas e carvao), com limites de P e Q e limite termico no CIA.


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
disp('.          AREA INTERCHANGE CONTROL: ACTIVSg2000, 6 AREAS          .')
disp('.                                                                  .')
disp('....................................................................')



%% ________________________________________________________________________
%%
%%                       2. CONFIGURACOES GLOBAIS
%% ________________________________________________________________________

%% a) Horas e Parques
horas  = [5368 8504];
nomeH  = {'Pico (11/08 15h)', 'Carga leve e eolica alta (20/12 07h)'};

% Numeracao original das barras
parque = [2127 4153];
metodos = {'HELM', 'NR'};



%% ________________________________________________________________________
%%
%%                       3. ESTUDO EM CADA HORA
%% ________________________________________________________________________

for ih = 1:numel(horas)

    %% a) Caso da Hora
    cfg      = Config_CIA('activsg_h');
    cfg.caso = sprintf('case_ACTIVSg2000_h%04d', horas(ih));
    mpc      = Carrega_Caso(cfg);
    wf       = zeros(1, numel(parque));

    for i = 1:numel(parque)
        wf(i) = find(mpc.numOrig == parque(i));
    end


    %% b) Ponto Base (sem CIA) e Intercambios Programados
    cfg0 = cfg;
    cfg0.NIesp_MW = zeros(1, 5);
    D0    = CIA_Dados(mpc, cfg0);
    tipo0 = D0.tipo;
    tipo0(D0.gci) = 2;
    V0    = CIA_Resolver(D0, 'NR', tipo0, D0.P0, D0.Q0, []);
    NI0   = D0.baseMVA*CIA_Intercambio(D0, V0);

    % Eolica prevista [MW]
    Peol  = D0.baseMVA*(D0.P0(wf) + D0.Pd(wf));
    [~, ~, cg0] = CIA_FluxoRamos(D0, V0);

    fprintf('\n===== %s =====\n', nomeH{ih});
    fprintf('Carga %.0f MW | eolica prevista: %.1f MW (barra %d) e %.1f MW (barra %d)\n', ...
        D0.baseMVA*sum(D0.Pd), Peol(1), parque(1), Peol(2), parque(2));
    fprintf('Base: V em [%.3f, %.3f] p.u. | carregamento maximo %.1f%%\n', ...
        min(abs(V0)), max(abs(V0)), 100*max(cg0));
    fprintf('Intercambios programados (recebidos) [MW]: %s\n', mat2str(NI0.', 6));


    %% c) Evento: Perda dos Parques
    cfg.NIesp_MW = NI0.';
    cfg.WF       = wf;
    cfg.WF_MW    = Peol.';

    % Injecao dos parques zerada
    D = CIA_Dados(mpc, cfg, [1, -1, -1]);


    %% d) Sem CIA
    tipoS = D.tipo;
    tipoS(D.gci) = 2;
    VS  = CIA_Resolver(D, 'NR', tipoS, D.P0, D.Q0, []);
    NIS = D.baseMVA*CIA_Intercambio(D, VS);
    [~, ~, cgS] = CIA_FluxoRamos(D, VS);
    [~, LS] = max(cgS);
    k   = D.ref;
    dPf = D.baseMVA*real(VS(k)*conj(D.Ybus(k,:)*VS) - V0(k)*conj(D.Ybus(k,:)*V0));

    fprintf('Sem CIA : desvio dos intercambios %s MW | folga +%.1f MW\n', ...
        mat2str((NIS - NI0).', 4), dPf);
    fprintf('          V em [%.3f, %.3f] p.u. | carregamento maximo %.1f%% (ramo %d-%d)\n', ...
        min(abs(VS)), max(abs(VS)), 100*cgS(LS), mpc.numOrig(D.ramo.de(LS)), mpc.numOrig(D.ramo.para(LS)));


    %% e) Com CIA (HELM e NR)
    for m = 1:2
        tic;
        r = CIA_Controle(D, metodos{m});
        t = toc;
        [~, ~, cg] = CIA_FluxoRamos(D, r.V);
        [~, L] = max(cg);

        fprintf('%-4s    : %d estagios, custo %d, %.1f s | desvio max %.1e MW\n', metodos{m}, ...
            r.estagios, r.custo, t, max(abs(D.baseMVA*r.NI - NI0)));
        fprintf('          lambda [MW] = %s | folga (Area 5) %.2f MW\n', mat2str(D.baseMVA*r.lam.', 5), D.baseMVA*r.lamRef);
        fprintf('          GCIs no limite de P: %d | barras com Q no limite: %d | ramos no limite: %d\n', ...
            numel(r.fixados), numel(r.fixQ), sum(r.termL > 0));
        fprintf('          V em [%.3f, %.3f] p.u. | carregamento maximo %.1f%% (ramo %d-%d)\n', ...
            min(abs(r.V)), max(abs(r.V)), 100*cg(L), mpc.numOrig(D.ramo.de(L)), mpc.numOrig(D.ramo.para(L)));

        R(ih).(metodos{m}) = r;
    end

    fprintf('max|V_HELM - V_NR| = %.2e p.u.\n', max(abs(R(ih).HELM.V - R(ih).NR.V)));

    R(ih).hora = horas(ih);
    R(ih).NI0  = NI0;
    R(ih).NIS  = NIS;
    R(ih).cgS  = max(cgS);
end

save('Resultados_ACTIVSg2000.mat', 'R');
