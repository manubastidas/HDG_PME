%%  RUN ALL
close all
clear
clc

%% INPUTS
mesh_list  = 1:3;
Orden_list = 1;
parm_list  = 2;

% Time parameters
T_MAX      = 1;
Initial_dt = 0.25;

n_mallas = length(mesh_list);
n_orden  = length(Orden_list);
n_parm   = length(parm_list);

dt_val  = [Initial_dt Initial_dt/4 Initial_dt/16 Initial_dt/64];
times_list = dt_val;

%% START PARALELLIZING

[ii,jj,kk] = meshgrid(parm_list,Orden_list,(mesh_list));
ii = permute(ii,[1 3 2]);
jj = permute(jj,[2 1 3]);
kk = permute(kk,[3  2 1]);
all_param =[ii(:) jj(:) kk(:)];

parm_list  = all_param(:,1);
Orden_list = all_param(:,2);
mesh_list  = all_param(:,3);
n_forloop  = size(all_param,1);

%% RUN !!!
for abc = 1:n_forloop
    
    % Select the parameters
    parm  = parm_list(abc,1);
    Orden = Orden_list(abc,1);
    mesh  = mesh_list(abc,1);
    
    Ncuad = 3*Orden;
    
    ArchivoL = sprintf('pre_mesh%iO%i',mesh,Orden);
    
    n_times      = times_list(mesh);
    TimetnSteps  = round(T_MAX/n_times);
    Timedt       = T_MAX./TimetnSteps;
    Timetime_vec = linspace(0,T_MAX,TimetnSteps+1);
      
    fprintf('\n m = %i | Orden = %i | Mesh = %i \n',parm,Orden,mesh)
    
    % Run the HDG code
    Main(ArchivoL,Orden,mesh,Timedt,TimetnSteps,parm)
end


