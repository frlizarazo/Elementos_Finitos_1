    %% defino las variables
    nef = 12;
    nno = 3*nef + 1;              % numero de nodos
    ngdl = nno;                   % el # de grados de libertad es el mismo # de nodos
    E   = 2e8;    % Pa          % modulo de elasticidad de la barra
    L   = 3;        % m           % longitud de la barra
    P = 5000;
    
    xnod = linspace(0,L,nno);     % posicion de los nodos
    le   = repmat(L/nef, nef, 1); % longitud de cada EF
    
    LaG=lagfunc(nef,4);

    %% Relacion de cargas puntuales
    f = zeros(nno,1); % vector de fuerzas nodales equivalentes global
    f(6*nef/L+1) = P;        % relaciono la carga puntual en el nodo "nno"

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
    rr1=zeros(1,nef);
    rr4=zeros(1,nef);
    for e = 1:nef      % ciclo sobre todos los elementos finitos
       idx  = LaG(e,:);
       r1   = r(3*e-2);
       r4   = r(3*e+1);
       rr1(e) = r1;
       rr4(e) =r4;
       [Ke,fe]=k_f_4n(nef,e,L,E,r1,r4);

       K(idx,idx) = K(idx,idx) + Ke;
       f(idx,:)   = f(idx,:)   + fe;
    end
    
    %% grados de libertad del desplazamiento conocidos y desconocidos
    c = [1 nno];    d = setdiff(1:nno,c);

    %% extraigo las submatrices y especifico las cantidades conocidas
    Kcd = K(c,d); fd = f(c); %Kcc = K(c,c); No se usan
    Kdd = K(d,d); fc = f(d); %Kdc = K(d,c);
    
    %% Areas
    A1=pi*rr1.^2;
    A4=pi*rr4.^2;

    %% resuelvo el sistema de ecuaciones
    ad = Kdd\fc;   % calculo desplazamientos desconocidos
    qd = Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
    a = zeros(ngdl,1);  a(d) = ad; % desplazamientos 
    q = zeros(ngdl,1);  q(c) = qd;             % fuerzas nodales equivalentes
    
    %% se realizan unos calculos intermedios que necesitaremos mas adelante
    nint = 10;           % numero de puntos donde se interpolará dentro del EF
    xi = linspace(-1,1,nint)'; % coordenadas naturales
    N = [- (9*xi.^3)/16+(9*xi.^2)/16+ xi/16-1/16,...
           (27*xi.^3)/16-(9*xi.^2)/16-(27*xi)/16+9/16,...
           (27*xi)/16-(9*xi.^2)/16-(27*xi.^3)/16+9/16,...
           (9*xi.^3)/16+(9*xi.^2)/16-xi/16-1/16];
    
    xx    = cell(nef,1); % interpol de posiciones (geometria) en el elemento
    uu    = cell(nef,1); % interpol desplazamientos en el elemento
    axial = cell(nef,1); % fuerzas axiales en el elemento
    esfuerzo = cell(nef,1);
    deformacion = cell(nef,1);
    
    for e = 1:nef       % ciclo sobre todas los elementos finitos
    
       Je = le(e)/2;     % Jacobiano del elemento ( = dx_dxi)
       Be = (1/Je)*[-(27*xi.^2)/16+(9*xi)/8+1/16,...
                     (81*xi.^2)/16-(9*xi)/8-27/16,...
                    -(81*xi.^2)/16-(9*xi)/8+27/16,...
                     (27*xi.^2)/16+(9*xi)/8-1/16]; % matriz de deformacion del elemento
    
       % vector de desplazamientos nodales del elemento a^{(e)}
       ae = [a(LaG(e,1));    a(LaG(e,2));    a(LaG(e,3)); a(LaG(e,4))]; % = a(LaG(e,:))';
    
       % vector de posiciones de nodos locales
       xe = [xnod(LaG(e,1)); xnod(LaG(e,2)); xnod(LaG(e,3)); xnod(LaG(e,4))]; %=xnod(LaG(e,:))';
    
       xx{e} = N*xe; % interpola sobre la geometría (coord naturales a geométricas)
       uu{e} = N*ae; % interpola sobre los desplazamientos
    
       N1_2n = (1-xi)/2;      N2_2n = (1+xi)/2;
       A = N1_2n.*A1 + N2_2n.*A4 ;
       Ae= A(:,e);
       De= E*Ae;
    
    %    Ae=[A(LaG(e,1)); A(LaG(e,2)); A(LaG(e,3)); A(LaG(e,4))];
    %    De=E*Ae;
       axial{e} = De.*Be*ae; % fuerza axial
       esfuerzo{e}=axial{e}./Ae; % esfuerzos
       deformacion{e}=(esfuerzo{e}./E)*100; %deformaciones
    end
    


    %% Grafico
    figure
    subplot(2,2,1)
    plot(xnod,a,'b')
    title('EF 4 nodos para los desplazamientos');
    xlabel('Eje X (m)')          % titulo del eje X
    ylabel('Desplazamientos (m)')    % titulo del eje Y
    legend('solucion por EF 4 nodos','Location','NorthEast');
    
    % 2) grafico los esfuerzos de la barra
    subplot(2,2,2);              % grafique en la parte inferior (2) del lienzo
    hold on;                     % no borre el lienzo
    for e = 1:nef % ciclo sobre todos los elementos finitos
        plot(xx{e}, axial{e}, 'r'); % grafico solucion por MEF
    end
    title('EF 4 nodos para la carga axial');
    xlabel('Eje X (m)')          % titulo del eje X
    ylabel('Carga axial (N)')    % titulo del eje Y
    legend('solucion por EF 4 nodos','Location','NorthEast');
    
    % 3) grafico los esfuerzos de la barra
    subplot(2,2,3);              % grafique en la parte inferior (2) del lienzo
    hold on;                     % no borre el lienzo
    for e = 1:nef % ciclo sobre todos los elementos finitos
        plot(xx{e}, esfuerzo{e}, 'g'); % grafico solucion por MEF
    end
    title('EF 4 nodos para los Esfuerzos');
    xlabel('Eje X (m)')          % titulo del eje X
    ylabel('Esfuerzos (N/m^2)')    % titulo del eje Y
    legend('solucion por EF 4 nodos','Location','NorthEast');
    
    % 4) grafico los esfuerzos de la barra
    subplot(2,2,4);              % grafique en la parte inferior (2) del lienzo
    hold on;                     % no borre el lienzo
    for e = 1:nef % ciclo sobre todos los elementos finitos
        plot(xx{e}, deformacion{e}, 'c'); % grafico solucion por MEF
    end
    title('EF 4 nodos para las Deformaciones');
    xlabel('Eje X (m)')          % titulo del eje X
    ylabel('Deformaciones (%)')    % titulo del eje Y
    legend('solucion por EF 4 nodos','Location','NorthEast');

    
    