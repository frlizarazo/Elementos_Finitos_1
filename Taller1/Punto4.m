%% Defino las variables
E = 2e8;    % Pa          % modulo de elasticidad de la barra
L = 3;        % m           % longitud de la barra    % N/m   

%% Solucion de la ecuacion diferencial

% Solucion numerica usando bvp4c (boundary value problem - MATLAB)
%   d /           du(x)  \
% ----| E(x) A(x)------- | + b(x) en x \in [0,L]     dado u(0) = 0
%  dx \            dx    /                                u(L) = 0
%

% En el caso mas general E, A y b son funciones. Escriba aqui las funciones
% como tal en caso de tener un caso mas general
EE = @(x) E;
AA = @(x) cono_Truncado(x);
bb = @(x) cargas(x);

% Se define la ecuacion diferencial, expresada como un sistema de dos
% ecuaciones diferenciales de primer orden
u      = 1; % y(1) = u(x)       
faxial = 2; % y(2) = faxial(x)
sist_eq_dif = @(x,y) [ y(faxial)/(EE(x)*AA(x)) 
                       -bb(x)                  ];

% Se definen las condiciones de frontera
% y_izq = condiciones de frontera del lado izquierdo (x=0)
% y_izq(1) = u(x=0)          y_izq(2) = faxial(x=0)
% y_der = condiciones de frontera del lado derecho   (x=L)
% y_der(1) = u(x=L)          y_der(2) = faxial(x=L)
cond_frontera = @ (y_izq,y_der) ...
                 [ y_izq(u)             % u(x=0) = 0 (desplazamiento)
                   y_der(u)];           % u(x=L) = 0 (carga axial)

% Solucion tentativa de la ecuacion diferencial
x = linspace(0,L,13);         % 30 puntos uniformemente distrib. entre 0 y L
y_inicial = bvpinit(x,[0 0]); % el [ 0 0 ] hace y_inicial.y = zeros(2,30)

% Solucion como tal de la ecuacion diferencial
sol = bvp4c(sist_eq_dif, cond_frontera, y_inicial);

% Evaluar la respuesta en los puntos x, ya que requiero evaluar la solucion
% diferentes puntos a los que retorna la funcion bvp4c()
y = deval(sol, x);


%% Grafico la solucion analitica y la solucion por el la funcion bvp4c
figure                       % cree un nuevo lienzo

% 1) grafico los desplazamientos de la barra
subplot(2,2,1);              % grafique en la parte superior (1) del lienzo
hold on;                     % no borre el lienzo 
plot(x, y(u,:), 'b');       % grafico solucion por bvp4c()
% plot(xx,a,'r');
title('La funcion bvp4c() para el desplazamiento');
xlabel('Eje X (m)')          % titulo del eje X
ylabel('Desplazamiento (m)') % titulo del eje Y
legend('solucion por bvp4c()', ...
                                                   'Location','SouthEast');

% 2) grafico la carga axial de la barra
subplot(2,2,2);              % grafique en la parte inferior (2) del lienzo
hold on;                     % no borre el lienzo
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
plot(x, y(faxial,:), 'r');  % grafico solucion por bvp4c()
title('La funcion bvp4c() para la carga axial');
xlabel('Eje X (m)')          % titulo del eje X
ylabel('Carga axial (N)')    % titulo del eje Y
legend('solucion por bvp4c()','Location','NorthEast');

% 3) grafico los esfuerzos de la barra
subplot(2,2,3);              % grafique en la parte inferior (2) del lienzo
hold on;                     % no borre el lienzo
esfuerzo=y(faxial,:);
for e=1:13
   esfuerzo(e)=esfuerzo(e)/cono_Truncado(x(e)); 
end
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
plot(x, esfuerzo, 'g');  % grafico solucion por bvp4c()
title('La funcion bvp4c() para la Esfuerzos');
xlabel('Eje X (m)')          % titulo del eje X
ylabel('Esfuerzos (N/m^2)')    % titulo del eje Y
legend('solucion por bvp4c()','Location','NorthEast');

% 4) grafico los esfuerzos de la barra
subplot(2,2,4);              % grafique en la parte inferior (2) del lienzo
hold on;                     % no borre el lienzo
deformacion=esfuerzo./E*100;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
plot(x, deformacion, 'c');  % grafico solucion por bvp4c()
title('La funcion bvp4c() para las deformaciones');
xlabel('Eje X (m)')          % titulo del eje X
ylabel('Deformaciones (%)')    % titulo del eje Y
legend('solucion por bvp4c()','Location','NorthEast');


