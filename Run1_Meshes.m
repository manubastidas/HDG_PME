%% PRE-PROCESS (COMPUTING MATRICES)
close all
clear
clc

% This code generates a set of regular meshes and pre-compute the local
% matrices necessary for the HDG implementation. 

%% INPUTS
mesh_list  = 1:4;
n_mallas   = length(mesh_list);

k_orden    = 1;
Orden_list = 1:k_orden;

div = [1/4 1/8 1/16 1/32 1/64];

for abc = 1:k_orden
    Orden = Orden_list(abc);
    Ncuad = 3*Orden;
    
    for pqr = 1:n_mallas
        
        % CREATE THE MESH
        mesh = mesh_list(pqr);
        [xx,yy]             = meshgrid(0:div(pqr):1,0:div(pqr):1);
        coordinate = [xx(:),yy(:)];
        element    = delaunay(xx,yy);
        
        nElement  = size(element,1);
        [cantFunciones,coefFunciones] = cardinalidad(Orden);
        
        [~,nodes2edge,~,edge2element,interioredge,GammaEdges]=edge(element,coordinate);
        
        hTriang = calcularh(element,coordinate,nElement);
        MaxH    = max(hTriang);
        disp(MaxH)
        
        % Compute de HDG matrices
        [ALocal,A,B,C,D,E,G,R,Afin,Jac,posF,posU,posV,posFlujo,edgesKnum,BASE] = ...
            SolucionHDG_prev(element,coordinate,Orden,cantFunciones,...
            coefFunciones,nodes2edge,edge2element,...
            interioredge,GammaEdges,Ncuad);
        
        % Save the pre-processing
        Archivo_save = sprintf('pre_mesh%iO%i',mesh,Orden);
        save(Archivo_save,'element','coordinate','nElement','cantFunciones','coefFunciones',...
            'nodes2edge','edge2element','interioredge','GammaEdges','MaxH','ALocal','A','B','C',...
            'D','E','G','R','Afin','Jac','posF','posU','posV','posFlujo','edgesKnum','BASE');
    end
end