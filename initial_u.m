function [u_0] = initial_u(nElement,coordinate,element,ALocal,Afin,BASE,parm)
% This function calculate the initial solution 

% Initialization
u_0 = zeros(BASE.cantFunciones*nElement,1);

M = ALocal(1:BASE.cantFunciones,1:BASE.cantFunciones);

for j=1:nElement
    pos   = BASE.cantFunciones*(j-1)+1:BASE.cantFunciones*j;
    coord = coordinate(element(j,:),:)';

    ZZ   = bsxfun(@plus,Afin{j,1}*([BASE.X1(:), BASE.Y1(:)]'),coord(:,1));
    u_ex = sol_exacta(ZZ(1,:)', ZZ(2,:)',0,parm);
    
    % Proyection of the initial solution
    Faux1 = repmat(u_ex,1,BASE.cantFunciones).*BASE.g_mn;
    for i= 1:BASE.cantFunciones
        F1(i,1) = BASE.WcuadX'*reshape(Faux1(:,i),size(BASE.WcuadX,1),...
            size(BASE.WcuadY,1))*BASE.WcuadY;
    end
   
    u_0(pos) = M\F1;
end