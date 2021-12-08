%% Solución Barra
%  Franklin Andrés Lizarazo Muñoz
%=========================================================================
% ||             b(x)             ||_ _             % E = 2e8 Pa
% ||\  ->->->->->->->->->->->->-> || |              % b(x) = 5x^2 kN/m
% || \____________________________|| |   _ _        %   
% ||                  -->5kN      || 9    | 3       % Sección Circular O
% ||  ____________________________|| cm  _|_cm      % 
% || /                            || |              % u(0) = 0
% ||/                             ||_|_             % u(L) = 0
% ||--------2m-------|            ||                %
% ||------------3m----------------||                %
%=========================================================================

%% Datos de entrada
nef  = 9;                     % numero de EF
nno  = 3*nef + 1;             % numero de nodos
ngdl = nno;                   % numero de grados de libertad
E    = 2e8;    % Pa           % modulo de elasticidad de la barra
L    = 3;      % m            % longitud de la barra
P    = 5000;   % N            % carga nodal al final de la barra

xnod = linspace(0,L,nno);     % posicion de los nodos
le   = repmat(L/nef, nef, 1); % longitud de cada EF

% Matriz de conectividad para EF de 4 nodos
LaG  = lagfunc(nef,4);

%% ensamblo la matriz de rigidez y el vector de f.n.e globales
syms xi x

xe = [0 1];               % Posiciones en las que cambia el radio
ri = [9/200 3/200];       % Radios dados
m  = diff(ri)/diff(xe);   % Pendiente de la función del radio
r  = zeros(nno,1);        % Vector de radios para cada nodo

% Relleno la medida de los radios en el vector
for n=1:nno
    if xnod(n)<1
        r(n)=r(n)+m*xnod(n)+ri(1);
    else
        r(n)=r(n)+ri(2);
    end
end

f = zeros(nno,1);        % vector de fuerzas nodales equivalentes global
f(6*nef/L+1) = P;        % relaciono la carga puntual ubicada a los 2 m

K  = zeros(ngdl);   % matriz de rigidez global

for e = 1:nef      % ciclo sobre todos los elementos finitos
   idx  = LaG(e,:);
   r1   = r(3*e-2);
   r4   = r(3*e+1);
   [Ke,fe]=k_f_4n(nef,e,L,E,r1,r4); %uso la función del punto 1

   K(idx,idx) = K(idx,idx) + Ke;   % ensamblo el K global
   f(idx,:)   = f(idx,:)   + fe;   % ensamblo el f global
end
%% grados de libertad restringidos y libres
c = [1 nno];    d = setdiff(1:nno,c);

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd |   | Kcc Kcd || ac |   | fd |
%|    | = |         ||    | - |    |
%| qc |   | Kdc Kdd || ad |   | fc |

%% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos
ac = 0;                 % desplazamientos conocidos

%% resuelvo el sistema de ecuaciones
ad = Kdd\fc;                               % calculo desplazamientos
qd = Kcd*ad - fd;                          % calculo fuerzas de equilibrio

a = zeros(ngdl,1);  a(c) = ac;  a(d) = ad; % desplazamientos 
q = zeros(ngdl,1);  q(c) = qd;             % fuerzas nodales equilibrio