function [B] = matrizB(cantFunciones,g_mn,g_mnX,g_mnY,Wx,Wy,Afin)
% This code calculate the submatrix B

B = zeros(2*cantFunciones,cantFunciones);
A = inv(Afin);

%% Integrals
for j = 1:cantFunciones
    for i=1:cantFunciones
        
        cadenaX = g_mnX(:,i)*A(1,1) + g_mnY(:,i)*A(2,1);
        cadenaY = g_mnX(:,i)*A(1,2) + g_mnY(:,i)*A(2,2);
        
        k = i+cantFunciones;
        
        % derX*funtion
        auxji = g_mn(:,j).*cadenaX;
        % derY*funtion
        auxjk = g_mn(:,j).*cadenaY;
        
        Terminoij = reshape(auxji,size(Wx,1),size(Wy,1));
        Terminokj = reshape(auxjk,size(Wx,1),size(Wy,1));
        
        % Integral !!!
        B(i,j)  = Wx'*Terminoij*Wy;
        B(k,j)  = Wx'*Terminokj*Wy;
    end
end