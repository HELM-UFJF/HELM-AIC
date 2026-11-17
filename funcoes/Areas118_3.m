function area = Areas118_3()
%AREAS118_3  Tres areas do IEEE 118 barras (zonas do arquivo de dados original).
%
%   Zona 1: barras 1-23, 25-32, 113-115, 117
%   Zona 2: barras 33-67
%   Zona 3: barras 24, 68-112, 116, 118 (contem a folga, barra 69)
%
%   Geradores reais (os demais sao compensadores sincronos, Pg = 0):
%     Zona 1: 10, 12, 25, 26, 31
%     Zona 2: 46, 49, 54, 59, 61, 65, 66
%     Zona 3: 69 (folga), 80, 87, 89, 100, 103, 111
%
%   Saida: area(k) = area da barra k (k = 1..118).


    %% a) Area de Cada Barra
    area = 3*ones(118, 1);
    area([1:23, 25:32, 113:115, 117]) = 1;
    area(33:67) = 2;
end
