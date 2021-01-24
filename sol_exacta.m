function u_ev = sol_exacta(xx,yy,tt,parm)
% This function calculates the exact solution to the PME 
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

