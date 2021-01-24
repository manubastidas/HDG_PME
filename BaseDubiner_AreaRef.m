function [g_mn,g_mnX,g_mnY,Wx1,Wy1,X1,Y1]= BaseDubiner_AreaRef(Ncuad,cantFunciones,coefFunciones)
% This code calculate the DUBINER BASIS that are use in the HDG
% implementation (2D)

%% Gauss - Quadrature on the reference triangle
% The order is 2N+1
[X1,Y1,Wx1,Wy1] = triquad(Ncuad,[0 0; 1 0; 0 1]); %N^2 points

% Some terms of the Jacobi polinomials
a1 = ((2*X1(:))./(1-Y1(:)))-1;
b1 = 2*Y1(:)-1;

%% Matrices
g_mn  = sparse(size(a1,1),cantFunciones);
g_mnX = sparse(size(a1,1),cantFunciones);
g_mnY = sparse(size(a1,1),cantFunciones);

for i = 1:cantFunciones
    
    % Index
    coef = coefFunciones(:,i);
    m = coef(1); n = coef(2);
    
    %% Dubiner Basis
    Terminouno0 = jacobiP(m,0,0,a1);
    Terminodos0 = jacobiP(n,2*m+1,0,b1);
    
    % Function
    g_mn(:,i) = 2^m.*Terminouno0.*(1-Y1(:)).^m.*Terminodos0;
    
    %% Gradient Dubiner basis
    if m == 0
        Terminouno1 = zeros(size(a1,1),1);
        Terminouno2 = zeros(size(a1,1),1);
    else
        Terminouno1 = jacobiP(m-1,1,1,a1);
        Terminouno2 = jacobiP(m-1,1,1,a1);
    end
    Terminodos1 = jacobiP(n,2*m+1,0,b1);
    Terminodos2  = jacobiP(n,2*m+1,0,b1);
    Terminotres2 = jacobiP(m,0,0,a1);
    if n == 0
        Terminocua2 = zeros(size(b1,1),1);
    else
        Terminocua2 = jacobiP(n-1,2*m+2,1,b1);
    end
    
    aux0 = 2^m.*(1-Y1(:)).^m.*((m+1)./(1-Y1(:)));
    % Gradient (component_x)
    g_mnX(:,i) = aux0.*Terminouno1.*Terminodos1;
    
    aux = (m+1)*(X1(:)./((1-Y1(:)).^2)).*Terminouno2;
    aux = aux.*(((1-Y1(:)).^m).*Terminodos2);
    aux2 = (-m*(1-Y1(:)).^(m-1)).*Terminodos2;
    aux2 = aux2 + ((1-Y1(:)).^m).*(n+2*m+2).*Terminocua2;
    aux2 = Terminotres2.*aux2;
    
    % Gradient (component_y)
    g_mnY(:,i) = 2^m*(aux+aux2);
end
