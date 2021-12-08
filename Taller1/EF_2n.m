function [q,a,f]=EF_2n(nef)
%Esto es una simplificacion de EF_4n, es en extremo similar
    %% defino las variables
    nno = nef + 1;              % numero de nodos
    ngdl = nno;                   % el # de grados de libertad es el mismo # de nodos
    E   = 2e8;    % Pa          % modulo de elasticidad de la barra
    L   = 3;        % m           % longitud de la barra
    P   = 5000;      % N           % carga nodal al final de la barra

    LaG=lagfunc(nef,2);
    
    %% Relacion de cargas puntuales
    f = zeros(nno,1); % vector de fuerzas nodales equivalentes global
    f(2/3*nef+1) = P;        % relaciono la carga puntual en el nodo "nno"

    %% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
    %  equivalentes global
    syms xi x
    xe=[0 1];      ri=[9/200 3/200];
    m=diff(ri)/diff(xe);
    r=zeros(nno,1);
    xx=linspace(0,3,nno);
    for n=1:nno
        if xx(n)<1
            r(n)=r(n)+m*xx(n)+ri(1);
        else
            r(n)=r(n)+ri(2);
        end
    end
    
    K  = zeros(ngdl);   % matriz de rigidez global
    for e = 1:nef      % ciclo sobre todos los elementos finitos
       idx  = LaG(e,:);
       r1   = r(e);
       
       [Ke,fe]=k_f_2n(nef,e,L,E,r1);

       K(idx,idx) = K(idx,idx) + Ke;
       f(idx,:)   = f(idx,:)   + fe;
    end

    %% grados de libertad del desplazamiento conocidos y desconocidos
    c = [1 nno];    d = setdiff(1:nno,c);

    %% extraigo las submatrices y especifico las cantidades conocidas
    Kcd = K(c,d); fd = f(c);
    Kdd = K(d,d); fc = f(d);

    %% resuelvo el sistema de ecuaciones
    ad = Kdd\fc;   % calculo desplazamientos desconocidos
    qd = Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
    a = zeros(ngdl,1);  a(d) = ad; % desplazamientos 
    q = zeros(ngdl,1);  q(c) = qd;             % fuerzas nodales equivalentes

    %% Grafico
    figure
    subplot(3,1,1)
    plot(xx,q,'r')
    
    subplot(3,1,2)
    plot(xx,a,'b')
    
    subplot(3,1,3)
    plot(xx,f,'g')
    
end