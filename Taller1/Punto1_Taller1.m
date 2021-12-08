% 1. Haga un programa para calcular la matriz de rigidez K y el vector de 
%    fuerzas nodales equivalentes f para elementos finitos de cuatro nodos
%    con sección transversal cónica truncada y bases de radio izquierdo r₁ 
%    (en x = x₁) y radio derecho r₂ (en x = x₂). Aquí se deben deducir las
%    fórmulas para este tipo especial de elemento finito.

%% Definicion de variables
syms xi x1 E L_e r1 r4

%% Defino las posiciones relativas de los 4 nodos
x4 = x1+3*L_e/3;                % Donde L_e es la longitud de cada EF
x3 = x1+2*L_e/3;                % y x1 seria la posición del ultimo nodo del
x2 = x1+L_e/3;                  % EF anterior (para el 1 EF x1=0)

%% Funciones de forma lagrangianas

% Se calculan los coeficientes de las funciones de forma
c1_4 = polyfit([-1 -1/3 1/3 1],[1 0 0 0],3);
c2_4 = polyfit([-1 -1/3 1/3 1],[0 1 0 0],3);
c3_4 = polyfit([-1 -1/3 1/3 1],[0 0 1 0],3);
c4_4 = polyfit([-1 -1/3 1/3 1],[0 0 0 1],3);

c1_3 = polyfit([-1 -1/3 1],[1 0 0],2);  
c2_3 = polyfit([-1 -1/3 1],[0 1 0],2);
c3_3 = polyfit([-1 -1/3 1],[0 0 1],2);

c1_2 = polyfit([-1 1],[1 0],1);
c2_2 = polyfit([-1 1],[0 1],1);

% Se eliminan los errores en la aproximacion numerica, haciendo los 
% coeficientes demasiado pequenios igual a cero
c1_4(abs(c1_4) < 1e-10) = 0;
c2_4(abs(c2_4) < 1e-10) = 0;
c3_4(abs(c3_4) < 1e-10) = 0;
c4_4(abs(c4_4) < 1e-10) = 0;

c1_3(abs(c1_3) < 1e-10) = 0;
c2_3(abs(c2_3) < 1e-10) = 0;
c3_3(abs(c3_3) < 1e-10) = 0;

% con los coeficientes corregidos se calculan las funciones de forma
N1_4 = poly2sym(c1_4, xi);    % -(9*xi^3)/16 + (9*xi^2)/16 + xi/16 - 1/16
N2_4 = poly2sym(c2_4, xi);    %  (27*xi^3)/16 - (9*xi^2)/16 - (27*xi)/16 + 9/16
N3_4 = poly2sym(c3_4, xi);    %  (27*xi)/16 - (9*xi^2)/16 - (27*xi^3)/16 + 9/16
N4_4 = poly2sym(c4_4, xi);    %  (9*xi^3)/8 + (9*xi^2)/8 - xi/8 - 1/8

N1_3 = poly2sym(c1_3, xi);    % Para interpolar la carga parabolica
N2_3 = poly2sym(c2_3, xi);
N3_3 = poly2sym(c3_3, xi);

N1_2 = poly2sym(c1_2, xi);    % Para interpolar el area
N2_2 = poly2sym(c2_2, xi);

%% Interpolacion de la geometria y sus derivadas
x   = simplify(N1_4*x1 + N2_4*x2 + N3_4*x3 + N4_4*x4);       
dx_dxi = diff(x, xi);

dxi_dx = 1/dx_dxi;
% recuerde que se debe garantizar que dx_dxi>0 y dxi_dx>0

%% Defino la carga distribuida
% Esta parte del codigo puede ser optimizada pero de momento para este
% taller la dejo asi

nef=3;
nno=nef*3+1;
L=3;
x=linspace(0,3,nno);

b=zeros(nno,1);     %Evaluo b para los nodos y almaceno en un vector
for n=1:nno
    b(n)=b(n)+5000*x(n)^2;
end

br=sym([]);
for n=1:nef
     b1=b(3*n-2);   %Defino 3 puntos sobre los cuales interpolar
     b2=b(3*n-1);
     b4=b(3*n+1);
     br(n)=N1_3*b1+N2_3*b2+N3_3*b4; %Interpolo con las funciones de forma
end

%% Definicion de la matriz de forma N y matriz de deformacion del elemento B
N = [N1_4 N2_4 N3_4 N4_4];
B = diff(N,xi)*dxi_dx;

%% Defino el Area en funcion de los radios 1 y 4
A1   = pi*r1^2;       A4   = pi*r4^2;
A    = N1_2*A1 + N2_2*A4;                  %Interpolo con las funciones de forma
    
%% Defino la matriz constitutiva
D = E*A;

%% Calculo la matriz de rigidez del elemento
K = int(B.'*D*B*dx_dxi, xi, -1, +1); K

%% Calculo la matriz de fuerzas nodales equivalentes del elemento
% En este caso como la carga parabolica va a tener diferente impacto sobre
% cada elemento finito hallo f para cada elemento finito, si bien se puede
% generalizar, como para el taller no lo vi necesario no lo hice.

b_e=br(1);
f1= int(N.'*b_e*dx_dxi, xi, -1, +1); f1

b_e=br(2);
f2= int(N.'*b_e*dx_dxi, xi, -1, +1); f2

b_e=br(3);
f3= int(N.'*b_e*dx_dxi, xi, -1, +1); f3

%% Defino f general para todos los EF
% por facilidad y paz mental este f lo obtuve en excel tras analizar las
% variaciones dentro de los f para cada EF sin embargo por el mismo sin
% sabor se haberlo sacado en excel, no lo incluyo en el desarrollo de los
% siguientes puntos, ademas de que no es necesario y lo pongo en esta parte
% unicamente por motivos esteticos para el informe.

syms n_nef nef
f_general = [5625*n_nef^2-9750*n_nef+4500
             16875*n_nef^2-27000*n_nef+10125
             16875*n_nef^2-6750*n_nef
             5625*n_nef^2-1500*n_nef+375].*L_e./nef^2;

% Aqui nef   = numero de elementos finitos que se estan usando y
%      n_nef = numero de elemento finito en el que se encuentra (1,2,...,nef)