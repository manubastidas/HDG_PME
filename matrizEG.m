function [E,G] = matrizEG(le,ne,g_mn_aristas,phi_arist,...
    cantFunciones,cantFunciones2D,W01,W0s2,Ncuad)
% This code calculate the submatrix E and G

% Change of the edges
le(1) = le(1)/sqrt(2);

% Inicialization
E = zeros(2*cantFunciones,3*cantFunciones2D);

% Sort the basis functions
punt = size(phi_arist,1);
aux_phi_arist = sparse(punt,3*cantFunciones2D);
for jj = 1:3
    aa  = (jj-1)*Ncuad+1; bb = jj*Ncuad;
    aa2 = (jj-1)*cantFunciones2D+1;
    bb2 =  jj*cantFunciones2D;
    aux_phi_arist(aa:bb,aa2:bb2) = phi_arist(aa:bb,:);
end
    
for i = 1:cantFunciones
    for k = 1:3*cantFunciones2D
        
        int = g_mn_aristas(:,i).*aux_phi_arist(:,k);
        int = reshape(int,Ncuad,3);
        
        integrale = sum(int.*[W0s2 W01 W01],1);
        E(i,k) = sum(le'.*integrale);
        E(i+cantFunciones,k) = E(i,k);
    end
end

% Auxiliar matrices with external normal
aux1_x = repmat(ne(1,:),cantFunciones2D,1); 
aux1_x = aux1_x(:);
aux1_x = aux1_x';
aux1_x = repmat(aux1_x,cantFunciones,1);

aux1_y = repmat(ne(2,:),cantFunciones2D,1);
aux1_y = aux1_y(:);
aux1_y = aux1_y';
aux1_y = repmat(aux1_y,cantFunciones,1);

G = E(1:cantFunciones,:);
E = E.*[aux1_x;aux1_y];

