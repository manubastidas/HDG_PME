function [R] = matrizR(le,phi_arist,W01,W0s2,cantFunciones2D,Ncuad)
% This code calculate the submatrix r

% Change variable of the edges
le(1) = le(1)/sqrt(2);

% Initialization
R = sparse(3*cantFunciones2D,3*cantFunciones2D);

% Sort PHI_ARIST
punt = size(phi_arist,1);
aux_phi_arist = sparse(punt,3*cantFunciones2D);

for jj = 1:3
    aa  = (jj-1)*Ncuad+1; bb = jj*Ncuad;
    aa2 = (jj-1)*cantFunciones2D+1;
    bb2 =  jj*cantFunciones2D;
    aux_phi_arist(aa:bb,aa2:bb2) = phi_arist(aa:bb,:);
end

%% TERMS 
% Outside of the diagonal (The basis is ortogonal but not on the edges 
for i=1:3*cantFunciones2D
    for j = 1:3*cantFunciones2D
        
        % Terms
        int = aux_phi_arist(:,i).*aux_phi_arist(:,j);
        int = reshape(int,Ncuad,3);
        
        % INTEGRAL!
        integrale = sum(int.*[W0s2 W01 W01],1);
        R(i,j) = sum(le'.*integrale);
    end
end