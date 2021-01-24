function [g_mn_aristas,w01,ws2,X2,Y2]= BaseDubiner_LineaRef(Ncuad,cantFunciones,coefFunciones)
% This code calculate the DUBINER BASIS that are use in the HDG
% implementation (1D)

%% Points over one edge
[x, w01]  = lgwt(Ncuad,0,1);
[x2,ws2]  = lgwt(Ncuad,0,sqrt(2));
t     = 3*pi/4;
vtra  = [cos(t) -sin(t);
    sin(t) cos(t)]*[x2';zeros(1,Ncuad)];

% figure
% scatter(x,zeros(2,1))
% hold on
% scatter(zeros(2,1),y)
% hold on
% scatter(vtra(1,:)'+1,vtra(2,:)')

%% Coordinates on the edges
X2 = [vtra(1,:)'+1 zeros(Ncuad,1) x];
Y2 = [vtra(2,:)' x zeros(Ncuad,1)];

X2(:,[1,3])=flip(X2(:,[1,3]));
Y2(:,[1,3])=flip(Y2(:,[1,3]));

a2 = ((2*X2(:))./(1-Y2(:)))-1;
b2 = 2*Y2(:)-1;

%% Matrices
g_mn_aristas  = zeros(size(a2,1),cantFunciones);

for i = 1:cantFunciones
    
    % Index
    coef = coefFunciones(:,i);
    m = coef(1); n = coef(2);
    
    %% Functions Dubiner basis 1D
    Terminouno0_ar = jacobiP(m,0,0,a2);
    Terminodos0_ar = jacobiP(n,2*m+1,0,b2);
    
    g_mn_aristas(:,i) = 2^m.*Terminouno0_ar.*(1-Y2(:)).^m.*Terminodos0_ar;
    
end