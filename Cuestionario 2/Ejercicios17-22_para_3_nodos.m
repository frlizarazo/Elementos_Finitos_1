% =========================================================================
% 
% Ejercicios 17 a 22 para 3 nodos Elementos Finitos 1
%
% Franklin Andrés Lizarazo Muñoz
% Estudiante Ing. Civil
%
% =========================================================================
% 
% Matriz de rigidez local y vector de fuerzas nodales 
% equivalentes local para EF de barra de 3 nodos
%
%        b(x)
% ||>->->->->->->->   >->->->->->->->
% ||*======*======*...*======*======*--> P
% ||<---- L  ---->|
%
% =========================================================================

% Se definene las variables
syms x1 xi L E A b

% Se calculan las funciones de forma
N1=poly2sym(round(polyfit([-1 0 1],[1 0 0],2),10),xi);
N2=poly2sym(round(polyfit([-1 0 1],[0 1 0],2),10),xi);
N3=poly2sym(round(polyfit([-1 0 1],[0 0 1],2),10),xi);
% Aqui la función round ayuda a eliminar errores de aproximacion
% es decir, elimina lo que es practicamente 0

% Se define la posicion de los nodos en funcion del primero
x2=x1+L/2;
x3=x1+L;

% Se crean los vectores de posiciones y funciones de forma
x=[x1; x2; x3];
N=[N1 N2 N3];

% Se interpola sobre la geometría
xe=N*x;

% Se calcula el Jacobiano
dN=diff(N,xi);
J=dN*x;

% Se calcula la matriz de deformaciones del elemento
B=dN/J;

% Se calcula la matriz constitutiva del elemento
D=E*A;

% Se calcula la matriz K local a partir del PTV
K=int(B'*D*B*J,xi,-1,1);
disp('K = (AE)/(3L)*')
pretty(K/((A*E)/(3*L')))

% Se calcula el vector f local a partir del PTV
f=int(N'*b*J,xi,-1,1);
disp('f = (Lb)/6*')
pretty(6*f/L/b)
