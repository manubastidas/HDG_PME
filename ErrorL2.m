function [ErrorLm1_Relativo,ErrorL2_Relativo] = ErrorL2(TimetnSteps,Timetime_vec,element,coordinate,Solucion_u,parm,BASE)
% This code calculate the Lm+1 and the L2 relative error 

nElement = size(element,1);

ErrorDif_m1   = zeros(nElement,TimetnSteps);
ErrorRef_p_m1 = zeros(nElement,TimetnSteps);
ErrorDif      = zeros(nElement,TimetnSteps);
ErrorRef_p    = zeros(nElement,TimetnSteps);

for kk=2:TimetnSteps+1
    for j=1:nElement
        pos   = BASE.cantFunciones*(j-1)+1:BASE.cantFunciones*j;
        coord = coordinate(element(j,:),:)';
        Afin  = coord(:,[2,3])-coord(:,[1,1]);
    
        ZZ=bsxfun(@plus,Afin*([BASE.X1(:), BASE.Y1(:)]'),coord(:,1));
        u_ex = sol_exacta(ZZ(1,:)', ZZ(2,:)',Timetime_vec(kk),parm);
        
        % Function to integrate
        Funcion  = abs(sum(repmat(Solucion_u(pos,kk)',size(BASE.g_mn,1),1).*BASE.g_mn,2) - u_ex);
        
        % Integral (difference) - LM+1 Error
        ErrorDif_m1(j,kk) = abs(det(Afin))*BASE.WcuadX'*...
            reshape(Funcion.^(parm+1),size(BASE.WcuadX,1),size(BASE.WcuadY,1))*BASE.WcuadY;
        % Integral (solution)- LM+1 Error
        ErrorRef_p_m1(j,kk) = abs(det(Afin))*BASE.WcuadX'*...
            reshape(abs(u_ex).^(parm+1),size(BASE.WcuadX,1),size(BASE.WcuadY,1))*BASE.WcuadY;
    
        % Integral (difference) - L2 Error
        ErrorDif(j,kk) = abs(det(Afin))*BASE.WcuadX'*...
            reshape(Funcion.^2,size(BASE.WcuadX,1),size(BASE.WcuadY,1))*BASE.WcuadY;
        % Integral (solution) - L2 Error
        ErrorRef_p(j,kk) = abs(det(Afin))*BASE.WcuadX'*...
            reshape(abs(u_ex).^2,size(BASE.WcuadX,1),size(BASE.WcuadY,1))*BASE.WcuadY;
    end
end

%% Relative Errors
ErrorLm1_Relativo = (sum(ErrorDif_m1,1)./sum(ErrorRef_p_m1,1));
ErrorL2_Relativo  = (sum(ErrorDif,1)./sum(ErrorRef_p,1));
