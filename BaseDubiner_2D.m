
function [phi_arist] = BaseDubiner_2D(Ncuad,cantFunciones2D)
% This code calculate the DUBINER BASIS that are use in the HDG
% implementation. This is the 2D Dubiner basis restricted to one edge

[x, ~]  = lgwt(Ncuad,0,1); %Points 0 1 
[x2,~]  = lgwt(Ncuad,0,sqrt(2)); % Points 0 S2

X_aux = [-sqrt(2)*x2+1; -2*x+1; -2*x+1];

phi_arist = zeros(length(X_aux),cantFunciones2D);
% Evaluation over one edge
for j = 1:cantFunciones2D
    % Base over the edges
    phi_arist(:,j) = jacobiP(j-1,0,0,X_aux);
end

