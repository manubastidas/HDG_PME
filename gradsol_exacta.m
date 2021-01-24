function [gradu_evx,gradu_evy] = gradsol_exacta(xx,yy,tt,parm)
% This function calculates the gradient of the exact solution of the PME 
% Barenbratt 

% d = 2;
x0 = [0.5,0.5];
Cb = 0.005;
t0 = 0.1;

t = tt+t0;
pnorm = @(x,y) sqrt((x-x0(1)).^2 + (y-x0(2)).^2);
part0 = t.^(1./(2*parm));
part1 = @(x,y) Cb- ((parm-1)./(4*parm^2))*(pnorm(x,y)./part0).^2;

beta = @(x,y)  t^(-1/parm).*max(part1(x,y),0).^(1/(parm-1));
u_ev = beta(xx,yy);

part2x = @(x,y) -1/(2*parm)*(1/(t^((parm+1)/parm)))*max(part1(x,y),0).^(1/(parm-1)).*(x-x0(1));
part2y = @(x,y) -1/(2*parm)*(1/(t^((parm+1)/parm)))*max(part1(x,y),0).^(1/(parm-1)).*(y-x0(2));

gradu_evx = part2x(xx,yy);
gradu_evy = part2y(xx,yy);

