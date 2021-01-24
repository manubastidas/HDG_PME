function [A] = matrizA(cantFunciones,g_mn,Wx,Wy)
% This code calculate the submatrix A

A = zeros(2*cantFunciones,2*cantFunciones);
for i = 1:cantFunciones
    % Terms to integrate
    Termino = reshape(g_mn(:,i).^2,size(Wx,1),size(Wy,1));
    % Integral !!!
    A(i,i)  = Wx'*Termino*Wy;
end
A(i+1:end,i+1:end) = A(1:i,1:i);
