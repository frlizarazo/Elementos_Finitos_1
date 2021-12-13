    %% Definicion de variables
    syms xi r1 r4
    nef=3; %Numero de elementos finitos
    n_nef=1; %Numero de elemento finito (e)
    L=3; %m
    E=2e8; %Pa
    L_e=L/nef;

    %% Defino las posiciones relativas de los 4 nodos
    x1 = 0+L_e*n_nef-1;
    x4 = x1+3*L_e/3;
    x3 = x1+2*L_e/3;
    x2 = x1+L_e/3;

    %% Funciones de forma lagrangianas

    % Se calculan los coeficientes de los polinomios
    c1_4 = polyfit([-1 -1/3 1/3 1],[1 0 0 0],3);
    c2_4= polyfit([-1 -1/3 1/3 1],[0 1 0 0],3);
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

    N1_3 = poly2sym(c1_3, xi);
    N2_3 = poly2sym(c2_3, xi);
    N3_3 = poly2sym(c3_3, xi);

    N1_2 = poly2sym(c1_2, xi);
    N2_2 = poly2sym(c2_2, xi);

%% Interpolacion de la geometria y sus derivadas
    x   = simplify(N1_4*x1 + N2_4*x2 + N3_4*x3 + N4_4*x4);       
    dx_dxi = diff(x, xi);
    dxi_dx = 1/dx_dxi;

    %% Carga Distribuida
    nno=nef*3+1;
    x=linspace(0,3,nno);
    
    b=zeros(nno,1);
    for n=1:nno
        b(n)=b(n)+5000*x(n)^2;
    end

    br=sym([]);
    
    for n=1:nef
        b1=b(3*n-2);
        b2=b(3*n-1);
        b4=b(3*n+1);
        br(n)=N1_3*b1+N2_3*b2+N3_3*b4;
    end
    b_e=br(n_nef);

    %% Definicion de la matriz de forma N y matriz de deformacion del elemento B
    N = [N1_4 N2_4 N3_4 N4_4];
    B = diff(N,xi)*dxi_dx;

    %% Area Variable
    A1   = pi*r1^2;
    A4   = pi*r4^2;
    A    = N1_2*A1 + N2_2*A4;
    
    %% "matriz constitutiva"
    D = E*A;

    %% Calculo la matriz de rigidez del elemento
    K = int(B.'*D*B*dx_dxi, xi, -1, +1); K

    %% Calculo la matriz de fuerzas nodales equivalentes del elemento
    f = int(N.'*b_e*dx_dxi, xi, -1, +1); f
