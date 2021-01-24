%% POST-PROCESS CODE
close all
clear
clc

%% INPUTS
mesh_list0  = 1:3;
Orden_list0 = 1;
parm_list0  = 2;

% Time
T_MAX      = 1;
Initial_dt = 0.25;

n_mallas = length(mesh_list0);
n_orden  = length(Orden_list0);
n_parm   = length(parm_list0);

dt_val  = [Initial_dt Initial_dt/4 Initial_dt/16 Initial_dt/64];
times_list = dt_val;

%% START PARALELLIZING

[ii,jj,kk]=meshgrid(parm_list0,Orden_list0,mesh_list0);
ii=permute(ii,[1 3 2]);
jj=permute(jj,[2 1 3]);
kk=permute(kk,[3  2 1]);
all_param =[ii(:) jj(:) kk(:)]

parm_list0  = all_param(:,1);
Orden_list0 = all_param(:,2);
mesh_list0  = all_param(:,3);
n_forloop  = size(all_param,1);

%% RUN !!
% Post-processing and the Errors 

for abc = 1:n_forloop
    
    parm  = parm_list0(abc,1);
    Orden = Orden_list0(abc,1);
    mesh  = mesh_list0(abc,1);
    
    Ncuad = 3*Orden;
    
    ArchivoL2 = sprintf('Sol_O%i_mesh%i_m%i',Orden,mesh,parm);
    load(ArchivoL2)
    
    Timetime_vec = linspace(0,T_MAX,TimetnSteps+1);
   
    % Post-processing
    [SolPlot_u,SolPlot_u_cont] = PostProcesamiento(TimetnSteps,coordinate,element,Solucion_u,BASE);
    fprintf('\n Postproces (Malla = %i de %i) \n',mesh,mesh_list0(end))
    % Error L
    [ErrorLm1_Relativo,ErrorL2_Relativo] = ErrorL2(TimetnSteps,Timetime_vec,...
        element,coordinate,Solucion_u,parm,BASE);
    fprintf('\n Error Lm1 %.3g \n',ErrorLm1_Relativo(end).^(1/(parm+1)))
    % Error H
    [ErrorH1_Relativo] = ErrorH1(TimetnSteps,Timetime_vec,Timedt,element,...
        coordinate,Solucion_F,parm,BASE);
    fprintf('\n Error H1 %.3g \n',ErrorL2_Relativo(end).^(1/2))
    % Negative mass
    [Nmass_total] = NMass(TimetnSteps,element,coordinate,Solucion_u,BASE);
    fprintf('\n Negative mass %.3g \n',Nmass_total(end))
    
    ArchivoL2 = sprintf('Sol_O%i_mesh%i_m%i',Orden,mesh,parm);
    save(ArchivoL2)
end

%% TABLE OF THE ERRORS 
for stu = 1:n_parm
    parm = parm_list0(stu);
    
    for abc = 1
        Tabla  = zeros(n_mallas,4);
        Orden = Orden_list0(abc);
        
        for pqr = 1:n_mallas
            
            mesh = mesh_list0(pqr);
            
            ArchivoL2 = sprintf('Sol_O%i_mesh%i_m%i',Orden,mesh,parm);
            load(ArchivoL2)
            
            % Information Table
            Tabla(pqr,1) = MaxH;
            Tabla(pqr,2) = Timedt;
            Tabla(pqr,3) = (ErrorLm1_Relativo(end)).^(1/(parm+1));
            Tabla(pqr,4) = (ErrorL2_Relativo(end)).^(1/2);
            Tabla(pqr,5) = sqrt(ErrorH1_Relativo(end));
            
            nEdgeInt     = size(interioredge,1);
            nEdgeGamma   = size(GammaEdges,1);
            nEdgeTotal   = nEdgeInt+nEdgeGamma;
            Tabla(pqr,6) = (Orden+1)*nEdgeTotal;
            Tabla(pqr,7) = mean(totalIter(2:end));
            
            Tabla(pqr,8) = Nmass_total(end);
        end
        
        ArchivoL2 = sprintf('Tab_O%i_m%i',Orden,parm);
        save(ArchivoL2)
    end
end

