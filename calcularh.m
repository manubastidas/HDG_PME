
function h = calcularh(element,coordinate,nElement)
% This function calculate the size of the mesh 
% This could be optimize but the mesh-structures of Matlab

h = zeros(nElement,1);
for j = 1:nElement
    % Coord (x;y) de cada uno de los vertices del tríangulo
    p = coordinate(element(j,:),:)';
    
    % Lenghts
    le1 = norm(p(:,1)-p(:,2));
    le2 = norm(p(:,2)-p(:,3));
    le3 = norm(p(:,3)-p(:,1));
    
    h(j) = max([le1,le2,le3]);
end