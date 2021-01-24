function [FF] = matrizF(Jac,ALocal,u_prev,u_it1,Llin,BASE,parm)
% This code compute the RHS for the HDG_PME 

suma = bsxfun(@times,BASE.g_mn,u_it1');

for i = 1:BASE.cantFunciones
    Aux = abs(suma(:,i)).^(parm-1).*u_it1(i).*BASE.g_mn(:,i).^2;
    % Terms to integrate
    Termino = reshape(Aux,size(BASE.WcuadX,1),size(BASE.WcuadY,1));
    % Integral !!!
    Faux11(i,1)  = BASE.WcuadX'*Termino*BASE.WcuadY;
end

% F1 and F2
Faux1    = u_prev.*diag(ALocal(1:BASE.cantFunciones,1:BASE.cantFunciones));
Faux12   = Llin.*u_it1.*diag(ALocal(1:BASE.cantFunciones,1:BASE.cantFunciones));

FF = [sparse(2*BASE.cantFunciones,1);Faux1*Jac;(Faux11-Faux12)*Jac];

end