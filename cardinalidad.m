function [cantFunciones,coef] = cardinalidad(Orden)
% This function calculate the size (cardinality of the DUBINER Basis)

% If zero order
if Orden == 0
    cantFunciones = 1;
    coef = [0 0]';
else
    % For order >= 1 
    coef  = CombVec(0:Orden,0:Orden);
    pos   = find(sum(coef,1)<=Orden);
    coef  = coef(:,pos);
    cantFunciones = length(pos);
end

end