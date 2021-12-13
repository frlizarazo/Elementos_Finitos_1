xx=linspace(0,3,13);
x=1:12;
a =[0.000000
    0.003884
    0.009661
    0.019048
    0.036749
    0.061801
    0.083376
    0.099955
    0.109741
    0.101820
    0.082686
    0.049713
    0.000000];

faxial =[16.562238
         16.471138
         16.145638
         15.429538
         14.166538
         12.200438
         9.374938
         5.533838
         -4.479162
         -10.820262
         -18.645762
         -28.111862];
     
Esfuerzos =[2.603e+003
            3.728e+003
            5.710e+003
            9.701e+003
            2.004e+004
            1.726e+004
            1.326e+004
            7.829e+003
            -6.337e+003
            -1.531e+004
            -2.638e+004
            -3.977e+004];
 
Deformaciones =Esfuerzos./E*100;

% 1) grafico los desplazamientos de la barra
subplot(2,2,1);              % grafique en la parte superior (1) del lienzo
hold on;                     % no borre el lienzo 
plot(xx,a, 'b');       % grafico solucion por Midas Gen
title(' Midas Gen para el desplazamiento');
xlabel('Eje X (m)')          % titulo del eje X
ylabel('Desplazamiento (m)') % titulo del eje Y
legend('solucion por Midas Gen', ...
                                                   'Location','SouthEast');

% 2) grafico la carga axial de la barra
subplot(2,2,2);              % grafique en la parte inferior (2) del lienzo
hold on;                     % no borre el lienzo
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
plot(x,faxial, 'r');  % grafico solucion por Midas Gen
title(' Midas Gen para la carga axial');
xlabel('Eje X (m)')          % titulo del eje X
ylabel('Carga axial (N)')    % titulo del eje Y
legend('solucion por Midas Gen','Location','NorthEast');

% 3) grafico los esfuerzos de la barra
subplot(2,2,3);              % grafique en la parte inferior (2) del lienzo
hold on;                     % no borre el lienzo
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
plot(x, Esfuerzos, 'g');  % grafico solucion por Midas Gen
title(' Midas Gen para la Esfuerzos');
xlabel('Eje X (m)')          % titulo del eje X
ylabel('Esfuerzos (N/m^2)')    % titulo del eje Y
legend('solucion por Midas Gen','Location','NorthEast');

% 4) grafico los esfuerzos de la barra
subplot(2,2,4);              % grafique en la parte inferior (2) del lienzo
hold on;                     % no borre el lienzo
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
plot(x, Deformaciones, 'c');  % grafico solucion por Midas Gen
title(' Midas Gen para las deformaciones');
xlabel('Eje X (m)')          % titulo del eje X
ylabel('Deformaciones (%)')    % titulo del eje Y
legend('solucion por Midas Gen','Location','NorthEast');

