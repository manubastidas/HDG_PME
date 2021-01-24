function [C] = matrizC(le,cantFunciones,g_mn_arista,W01,W0s2,Ncuad)
% This code calculate the submatrix C

% Change variable of the edges
le(1) = le(1)/sqrt(2); 
C     = zeros(cantFunciones,cantFunciones);

for i = 1:cantFunciones
    for j = i:cantFunciones
        
        % Terms to integrate (all points)
        int = g_mn_arista(:,i).*g_mn_arista(:,j);
        int = reshape(int,Ncuad,3);
        
        %% Integral!!!!
        integrale = sum(int.*[W0s2 W01 W01],1);

        C(i,j) = sum(le'.*integrale);
        C(j,i) = C(i,j);
    end
end

