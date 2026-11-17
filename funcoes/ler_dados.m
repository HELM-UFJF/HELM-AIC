function [De, Para, R, X, Z, Y, b, SHL, a, SHB, Tipo, V, V_esp, Teta, Qmin, Qmax, Pg, Pd, Qg, Qd, P_esp, Q_esp, Nbarra, Nlinha] = ler_dados(baseMVA, bus, gen, branch, FE, FC)
%LER_DADOS  Le os dados de barras e ramos de um caso no formato do codigo.
%
%   Entradas: baseMVA, bus, gen, branch e os fatores FE (injecoes) e FC
%             (cargas)
%   Saidas  : dados dos ramos (De, Para, R, X, Z, Y, b, SHL, a), das barras
%             (SHB, Tipo, V, V_esp, Teta, Qmin, Qmax, Pg, Pd, Qg, Qd),
%             injecoes especificadas (P_esp, Q_esp) e dimensoes do sistema.
%             Grandezas em p.u. e angulos em rad.



    %% a) Dados dos Ramos
    % Barras de origem e de destino
    De   = branch(:,1);
    Para = branch(:,2);

    % Resistencia, reatancia, impedancia e admitancia serie [p.u.]
    R = branch(:,3);
    X = branch(:,4);
    Z = R + 1j*X;
    Y = 1./Z;

    % Metade da susceptancia de carregamento e shunt correspondente
    b   = branch(:,5)/2;
    SHL = 1j*b;

    % Tap dos transformadores
    a = branch(:,9);


    %% b) Dados das Barras
    % Shunt da barra (bancos de capacitores ou reatores)
    SHB = 1j*bus(:,6)/baseMVA;

    % Tipo: 1 = PQ, 2 = PV, 3 = referencia
    Tipo = bus(:,2);

    % Tensao inicial e especificada [p.u.] e angulo inicial [rad]
    V     = bus(:,8);
    V_esp = bus(:,8);
    Teta  = bus(:,9)*pi/180;

    % Limites de potencia reativa [p.u.]
    Qmin = bus(:,16)/baseMVA;
    Qmax = bus(:,17)/baseMVA;

    % Geracao e carga [p.u.]
    Pg = bus(:,14)/baseMVA;
    Pd = FC.*bus(:,3)/baseMVA;
    Qg = bus(:,15)/baseMVA;
    Qd = FC.*bus(:,4)/baseMVA;


    %% c) Injecoes Especificadas
    P_esp = FE.*(Pg - Pd);
    Q_esp = FE.*(Qg - Qd);


    %% d) Dimensoes do Sistema
    Nbarra = length(bus(:,1));
    Nlinha = length(branch(:,1));
end
