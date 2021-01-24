function [Nmass_total] = NMass(TimetnSteps,element,coordinate,Solucion_u,BASE)
% This code calculate the negative mass of the solution 

nElement = size(element,1);

Int_up   = zeros(nElement,TimetnSteps);
Int_down = zeros(nElement,TimetnSteps);

for kk=1:TimetnSteps+1
    for j=1:nElement
        pos   = BASE.cantFunciones*(j-1)+1:BASE.cantFunciones*j;
        coord = coordinate(element(j,:),:)';
        Afin  = coord(:,[2,3])-coord(:,[1,1]);
    
        % Function
        Funcion  = sum(repmat(Solucion_u(pos,kk)',size(BASE.g_mn,1),1).*BASE.g_mn,2);
        Funcion_men = max(-Funcion,0);
        
        % Integral (Function)
        Int_up(j,kk)  = abs(det(Afin))*BASE.WcuadX'*reshape(Funcion_men,size(BASE.WcuadX,1),size(BASE.WcuadY,1))*BASE.WcuadY;
        % Integral -(Function)
        Int_down(j,kk)= abs(det(Afin))*BASE.WcuadX'*reshape(abs(Funcion),size(BASE.WcuadX,1),size(BASE.WcuadY,1))*BASE.WcuadY;
     end
end
%% Relative negative mass
Nmass_total = (sum(Int_up,1)./sum(Int_down,1));


