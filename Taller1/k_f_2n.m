function [K,f]=k_f_2n(nef,n_nef,L,E,r1,r4)
% Esto es una simplificación de la funcion k_f_4n.m

    %% Definicion de variables
    syms xi
    L_e=L/nef;

    %% Defino las posiciones relativas de los 4 nodos
    x1 = 0+L_e*n_nef-1;
    x2 = x1+L_e;

    %% Funciones de forma lagrangianas

    % Se calculan los coeficientes de los polinomios
    c1_3 = polyfit([-1 -1/3 1],[1 0 0],2);
    c2_3 = polyfit([-1 -1/3 1],[0 1 0],2);
    c3_3 = polyfit([-1 -1/3 1],[0 0 1],2);

    c1_2 = polyfit([-1 1],[1 0],1);
    c2_2 = polyfit([-1 1],[0 1],1);

% Se eliminan los errores en la aproximacion numerica, haciendo los 
% coeficientes demasiado pequenios igual a cero
    c1_3(abs(c1_3) < 1e-10) = 0;
    c2_3(abs(c2_3) < 1e-10) = 0;
    c3_3(abs(c3_3) < 1e-10) = 0;

    % con los coeficientes corregidos se calculan las funciones de forma
    N1_3 = poly2sym(c1_3, xi);
    N2_3 = poly2sym(c2_3, xi);
    N3_3 = poly2sym(c3_3, xi);

    N1_2 = poly2sym(c1_2, xi);
    N2_2 = poly2sym(c2_2, xi);

%% Interpolacion de la geometria y sus derivadas
    x   = simplify(N1_2*x1 + N2_2*x2);       
    dx_dxi = diff(x, xi);
    dxi_dx = 1/dx_dxi;

    %% Carga Distribuida
    nno=3*nef+1;            %En esta parte se toman 4 nodos por EF con el
    x=linspace(0,3,nno);    %proposito de hacer una interpolación parabolica
                            %de la carga distribuida.
    b=zeros(nno,1);
    for n=1:nno
        b(n)=b(n)+5000*x(n)^2;
    end

    br=sym([]); %No es la mejor practica porque el vector se actualiza de tamaño en cada iteración.
    for n=1:nef
        b1=b(3*n-2);        %como se ve aca solo son necesarios 3 puntos
        b2=b(3*n-1);        %pero por la formulación del numeral 1 era 
        b4=b(3*n+1);        %conveniente dividir en 4 nodos.
        br(n)=N1_3*b1+N2_3*b2+N3_3*b4;
    end
    b_e=br(n_nef);

    %% Definicion de la matriz de forma N y matriz de deformacion del elemento B
    N = [N1_2 N2_2];
    B = diff(N,xi)*dxi_dx;

    %% Area Variable
    A   = pi*r1^2;
    
    %% "matriz constitutiva"
    D = E*A;

    %% Calculo la matriz de rigidez del elemento
    K = int(B.'*D*B*dx_dxi, xi, -1, +1);

    %% Calculo la matriz de fuerzas nodales equivalentes del elemento
    f = int(N.'*b_e*dx_dxi, xi, -1, +1);

end