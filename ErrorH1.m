function [ErrorH1_Relativo] = ErrorH1(TimetnSteps,Timetime_vec,Timedt,element,coordinate,Solucion_F,parm,BASE)
% This code calculate the H1 relative error 

nElement   = size(element,1);
Solucion_F = -Solucion_F/Timedt;

for kk=2:TimetnSteps+1
    parfor j=1:nElement
        posdux = 2*BASE.cantFunciones*(j-1)+1:BASE.cantFunciones*(2*j-1);
        posduy = BASE.cantFunciones*(2*j-1)+1:2*BASE.cantFunciones*j;
        
        coord = coordinate(element(j,:),:)';
        Afin  = coord(:,[2,3])-coord(:,[1,1]);
        
        % Transform to a basis functions
        ZZ=bsxfun(@plus,Afin*([BASE.X1(:), BASE.Y1(:)]'),coord(:,1));
        % Calculate gradient
        [gradu_evx,gradu_evy] = gradsol_exacta(ZZ(1,:)', ZZ(2,:)',Timetime_vec(kk),parm);
        
        % Function to integrate
        FuncionX  = abs(sum(repmat(Solucion_F(posdux,kk)',size(BASE.g_mn,1),1).*BASE.g_mn,2)...
            - gradu_evx);
        FuncionY  = abs(sum(repmat(Solucion_F(posduy,kk)',size(BASE.g_mn,1),1).*BASE.g_mn,2)...
            - gradu_evy);
        
        % Integrate (difference)
        ErrorDif(j,kk)=abs(det(Afin))*BASE.WcuadX'*...
            reshape(FuncionX.^2+FuncionY.^2,size(BASE.WcuadX,1),size(BASE.WcuadY,1))*BASE.WcuadY;
        
        % Integrate (Function)
        ErrorRef_g(j,kk)=abs(det(Afin))*BASE.WcuadX'...
            *reshape(gradu_evx.^2 + gradu_evy.^2 ,...
            size(BASE.WcuadX,1),size(BASE.WcuadY,1))*BASE.WcuadY;
    end
end
%% Relative error
ErrorH1_Relativo = (sum(ErrorDif,1)./sum(ErrorRef_g,1));