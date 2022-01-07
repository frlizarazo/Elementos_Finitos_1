% =========================================================================
% 
% Ejercicios 17 a 22 Elementos Finitos 1
%
% Franklin Andrés Lizarazo Muñoz
% Estudiante Ing. Civil
%
% =========================================================================
% 
% Matriz de rigidez local y vector de fuerzas nodales 
% equivalentes local para EF de barra de 4 nodos
%
%             b(x)
% ||->->->->->->->->->->->   ->->->->->->->->->->->
% ||*======*======*======*...*======*======*======*--> P
% ||<-------- L -------->|
%
% =========================================================================

% Se definene las variables
syms x1 xi L E A b

% Se calculan las funciones de forma
N1 = poly2sym(polyfit([-1 -1/3 1/3 1],[1 0 0 0],3),xi);
N2 = poly2sym(polyfit([-1 -1/3 1/3 1],[0 1 0 0],3),xi);
N3 = poly2sym(polyfit([-1 -1/3 1/3 1],[0 0 1 0],3),xi);
N4 = poly2sym(polyfit([-1 -1/3 1/3 1],[0 0 0 1],3),xi);

% Se define la posicion de los nodos en funcion del primero
x2 = x1+L/3;
x3 = x1+2/3*L;
x4 = x1+L;

% Se crean los vectores de posiciones y funciones de forma
x = [x1; x2; x3; x4];
N = [N1 N2 N3 N4];

% Se interpola sobre la geometría
xe = N*x;

% Se calcula el Jacobiano
dN = diff(N,xi);
J  = dN*x;

% Se calcula la matriz de deformaciones del elemento
B = dN/J;

% Se calcula la matriz constitutiva del elemento
D = E*A;

% Se calcula la matriz K local a partir del PTV
K = int(B'*D*B*J,xi,-1,1);
disp('K = (AE)/(40L)*')
pretty(K/((A*E)/(40*L')))

% Se calcula el vector f local a partir del PTV
f = int(N'*b*J,xi,-1,1);
disp('f = (Lb)/8*')
pretty(f/(L*b/8))