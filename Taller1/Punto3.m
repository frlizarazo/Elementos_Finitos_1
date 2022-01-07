    %% defino las variables
    nef  = 12;
    nno  = nef + 1;              % numero de nodos
    ngdl = nno;                   % el # de grados de libertad es el mismo # de nodos
    E    = 2e8;    % Pa          % modulo de elasticidad de la barra
    L    = 3;        % m           % longitud de la barra

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
    q = zeros(ngdl,1);  q(c) = qd; % fuerzas nodales equivalentes
    
    %% fuerza axial
    faxial = zeros(nef,1);
    for e = 1:nef % ciclo sobre todas los elementos finitos
        Be = [-nef/L nef/L];
        ae = a(LaG(e,:));
        A  = pi*r(e)^2;
        faxial(e) = (E*A)*Be*ae; % = D*B(e)*a(e)
    end
    
    %% Esfuerzos
    esfuerzo=zeros(nef,1);
    for e = 1:nef % ciclo sobre todas los elementos finitos
        Be = [-nef/L nef/L];
        ae = a(LaG(e,:));
        esfuerzo(e) = E*Be*ae; % = E*B(e)*a(e)
    end
    
    %% Deformaciones
    deformacion=esfuerzo./E*100;
    
    %% Grafico
    figure
    subplot(2,2,1)
    plot(xx,a,'b')
    title('EF 2 nodos para los desplazamientos');
    xlabel('Eje X (m)')          % titulo del eje X
    ylabel('Desplazamientos (m)')    % titulo del eje Y
    legend('solucion por EF 2 nodos','Location','NorthEast');
    
    subplot(2,2,2)
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    for e = 1:nef % ciclo sobre todas los elementos finitos
        hold on
        plot([xx(e) xx(e+1)], [faxial(e) faxial(e)], 'r.-'); % grafico solucion por MEF
    end
    title('EF 2 nodos para la fuerza axial');
    xlabel('Eje X (m)')          % titulo del eje X
    ylabel('Fuerza Axial (N)')    % titulo del eje Y
    legend('solucion por EF 2 nodos','Location','NorthEast');
    
    subplot(2,2,3)
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    for e = 1:nef % ciclo sobre todas los elementos finitos
        hold on
        plot([xx(e) xx(e+1)], [esfuerzo(e) esfuerzo(e)], 'g.-'); % grafico solucion por MEF
    end
    title('EF 2 nodos para laos esfuerzos');
    xlabel('Eje X (m)')          % titulo del eje X
    ylabel('Esfuerzos (N/m^2)')    % titulo del eje Y
    legend('solucion por EF 2 nodos','Location','NorthEast');
    
    subplot(2,2,4)
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    for e = 1:nef % ciclo sobre todas los elementos finitos
      hold on
      plot([xx(e) xx(e+1)], [deformacion(e) deformacion(e)], 'c.-'); % grafico solucion por MEF
    end
    title('EF 2 nodos para las deformaciones');
    xlabel('Eje X (m)')          % titulo del eje X
    ylabel('Deformaciones (%)')    % titulo del eje Y
    legend('solucion por EF 2 nodos','Location','NorthEast');


    