function [ALocal,A,B,C,D,E,G,R,Afin,Jac,posF,posU,posV,posFlujo,edgesKnum,BASE] = ...
    SolucionHDG_prev(element,coordinate,Orden,cantFunciones,coefFunciones,nodes2edge,edge2element,...
    interioredge,GammaEdges,Ncuad)

% Assemble matrices A - B - C
nElement   = size(element,1);
nEdgeInt   = size(interioredge,1);
nEdgeGamma = size(GammaEdges,1);
nEdgeTotal = nEdgeInt+nEdgeGamma;

A  = cell(nElement,1);
B  = cell(nElement,1);
C  = cell(nElement,1);
D = cell(nElement,1);
Afin = cell(nElement,1);
Jac  = zeros(nElement,1);

E = cell(nElement,1);
G = cell(nElement,1);
R = cell(nElement,1);

%% Index of the local positions
posF   = cell(nElement,1);
posU  = cell(nElement,1);
posV  = cell(nElement,1);

posFlujo  = cell(nElement,1);
edgesKnum = cell(nElement,1);

%%
BASE = struct();
BASE.cantFunciones = cantFunciones;
BASE.coefFunciones = coefFunciones;
[BASE.g_mn,BASE.g_mnX,BASE.g_mnY,BASE.WcuadX,...
    BASE.WcuadY,BASE.X1,BASE.Y1]=...
    BaseDubiner_AreaRef(Ncuad,BASE.cantFunciones,BASE.coefFunciones);

[BASE.g_mn_aristas,BASE.W01,BASE.W0s2,BASE.X2,BASE.Y2]= ...
    BaseDubiner_LineaRef(Ncuad,BASE.cantFunciones,BASE.coefFunciones);

BASE.cantFunciones2D = Orden +1;
[BASE.phi_arist] = BaseDubiner_2D(Ncuad,BASE.cantFunciones2D);

% Local submatrix A on the reference triangle
[ALocal] = matrizA(cantFunciones,BASE.g_mn,BASE.WcuadX,BASE.WcuadY);

parfor elemento = 1:nElement
    
    % Local positions (to ensamble)
    posF{elemento,1} = 2*cantFunciones*elemento-2*cantFunciones+1:2*cantFunciones*elemento;
    posV{elemento,1} = cantFunciones*elemento-cantFunciones+1:cantFunciones*elemento;
    posU{elemento,1} = cantFunciones*elemento-cantFunciones+1:cantFunciones*elemento;
    
    % Coord (x;y) - Each vertex
    coord = coordinate(element(elemento,:),:)';
    Afin{elemento,1}  = coord(:,[2,3])-coord(:,[1,1]);
    Jac(elemento,1) = abs(det(Afin{elemento,1}));
    
    %% MATRIX A
    A{elemento,1}= Jac(elemento,1)*ALocal;
    
    %% MATRIX B
    [BLocal]      = matrizB(cantFunciones,BASE.g_mn,BASE.g_mnX,BASE.g_mnY,...
        BASE.WcuadX,BASE.WcuadY,Afin{elemento,1});
    B{elemento,1} = Jac(elemento,1)*BLocal;
    
    %% MATRIX C
    edgesk = diag(nodes2edge(element(elemento,[2 3 1]),element(elemento,[3 1 2])));
    p  = [coord(:,[2 3 1]);coord(:,[3 1 2])];
    % Lenght all the edges
    le = [norm(p(1:2,1)-p(3:4,1));norm(p(1:2,2)-p(3:4,2));norm(p(1:2,3)-p(3:4,3))];
    
    C{elemento,1}  = matrizC(le,cantFunciones,BASE.g_mn_aristas,BASE.W01,BASE.W0s2,Ncuad);
    
    %% Matrix D
    D{elemento,1} = Jac(elemento,1)*ALocal(1:cantFunciones,1:cantFunciones);
    
    %% MATRIX E and G
    
    % Normal vector ->> 1:2 nx:ny
    normalvector = zeros(2,3);
    for n=1:3
        normalvector(1:2,n) = (p(1:2,n)-p(3:4,n))'*[0,-1;1,0]/le(n);
    end
    % Directions of edges on the triangle
    posicionesElem = zeros(3,BASE.cantFunciones2D);
    for ii=1:3
        trin = BASE.cantFunciones2D*edgesk(ii)-BASE.cantFunciones2D+1;
        posicionesElem(ii,:) = trin:BASE.cantFunciones2D*edgesk(ii);
    end
    posFlujo{elemento,1}  = posicionesElem';
    %     posFlujo{elemento,1}  = [2*edgesk'-1;2*edgesk'];
    %     OrdenposFlujo         = posFlujo{elemento,1}(:);
    
    % Local enumeration of the edges
    edgesKnum{elemento,1} = edgesk;
    
    % Orientation of the edges 
    aux = zeros(0);
    aux(element(elemento,:))=1:3;
    aaa = edge2element(edgesk,1:2);
    % (local) InitialPoint -> FinalPoint
    aaa= aux(aaa);
    % Standar orientation: 2 3 1
    % Cambios -> 1 If some edge has different orientation
    cambios = find(([2 3 1]'-aaa(:,1))~=0);
    
    % Reorganize the basis
    g_mn_aristas_aux=BASE.g_mn_aristas;
    for ii = 1:length(cambios)
        postochange = (cambios(ii)-1)*Ncuad+1:cambios(ii)*Ncuad;
        g_mn_aristas_aux(postochange,:) = flip(g_mn_aristas_aux(postochange,:));
    end
    
    [E{elemento,1}, G{elemento,1}] = matrizEG(le,normalvector,g_mn_aristas_aux,...
        BASE.phi_arist,cantFunciones,BASE.cantFunciones2D,BASE.W01,BASE.W0s2,Ncuad);
    E{elemento,1} = -E{elemento,1};
    
    %% MATRIX R
    R{elemento,1}= matrizR(le,BASE.phi_arist,BASE.W01,BASE.W0s2,...
        BASE.cantFunciones2D,Ncuad); 
    
end
