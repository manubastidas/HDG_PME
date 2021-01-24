function [SolPlot_u,SolPlot_u_cont] = PostProcesamiento(TimetnSteps,coordinate,element,Solucion_u,BASE)
% This code calculate the Postprocessing of the solutio  of the solution u

%% DUBINER BASIS -> Baricenter
bari = [0 0; 1  0; 0 1];
X1 = bari(:,1); Y1= bari(:,2);

% Termins Jacobi
a1 = (2*X1(:))./((1-Y1(:)))-1;
a1(isnan(a1))=-1;
b1 = 2*Y1(:)-1;
% Inicialization
g_mnBari  = sparse(size(a1,1),BASE.cantFunciones);
for i = 1:BASE.cantFunciones
    % Index
    coef = BASE.coefFunciones(:,i);
    m = coef(1); n = coef(2);
    
    % Dubiner terms
    Terminouno0 = jacobiP(m,0,0,a1);
    Terminodos0 = jacobiP(n,2*m+1,0,b1);
    %Function
    g_mnBari(:,i) = 2^m.*Terminouno0.*(1-Y1(:)).^m.*Terminodos0;
end

%% POST-PROCESS PRESSURE
nElement = size(element,1);
contador       = zeros(size(coordinate,1),1);
SolPlot_u_cont = zeros(size(coordinate,1),TimetnSteps+1);

for kk=1:TimetnSteps+1
    for j = 1:nElement
        pos  = BASE.cantFunciones*j-BASE.cantFunciones+1:BASE.cantFunciones*j;
        pos2 = 3*j-2:3*j;
        
        % Solution to plot
        SolPlot_u(pos2,kk) = g_mnBari*Solucion_u(pos,kk);
        
        SolPlot_u_cont(element(j,:),kk) = SolPlot_u_cont(element(j,:),kk)+SolPlot_u(pos2,kk);
        contador(element(j,:))       = contador(element(j,:)) + ones(3,1);
    end
    SolPlot_u_cont(:,kk) = SolPlot_u_cont(:,kk)./contador;
end



