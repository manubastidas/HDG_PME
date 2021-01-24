% This script make the plots of the convergence 
close all

load('Tab_O2_m4.mat')

OrdenErrorLm1= abs(log(Tabla(1:(end-1),3)./Tabla(2:(end),3)))'./abs(log(Tabla(1:(end-1),1)./Tabla(2:(end),1)))';
OrdenErrorL2 = abs(log(Tabla(1:(end-1),4)./Tabla(2:(end),4)))'./abs(log(Tabla(1:(end-1),1)./Tabla(2:(end),1)))';
OrdenErrorH1 = abs(log(Tabla(1:(end-1),5)./Tabla(2:(end),5)))'./abs(log(Tabla(1:(end-1),1)./Tabla(2:(end),1)))';

convLm1      = ceil(max(OrdenErrorLm1)*10)/10;
reference1 =  (Tabla(:,1)).^(convLm1);
convL2      = ceil(max(OrdenErrorL2)*10)/10;
reference2 =  (Tabla(:,1)).^(convL2);
convLH1      = ceil(max(OrdenErrorH1(end))*10)/10;
reference3 =  (Tabla(:,1).^(convLH1));

%% CONVERGENCE LM1
figure
hold on
plot(Tabla(:,1),Tabla(:,3),'-ob','linewidth',2,'MarkerSize',12)
plot(Tabla(:,1),Tabla(:,4),'-or','linewidth',2,'MarkerSize',12)
plot(Tabla(:,1),reference1,'-.k','linewidth',2)
plot(Tabla(:,1),reference2,'--k','linewidth',2)
set(gca, 'YScale', 'log','XScale', 'log')
set(gca,'fontsize',14)
set(gcf, 'Color', 'w')
box on
axis tight

ax = gca;
ax.YAxis.Exponent = 2;
% NumTicks = 4;
% L = get(gca,'YLim');
% set(gca,'YTick',linspace(ceil(L(1)*10)/10,ceil(L(2)*10)/10,NumTicks))
xtickformat('%.2f')
% ytickformat('%g')
xtickangle(45)

grid on
xlabel('h','fontsize',20,'Interpreter','latex')
ylabel('Error','fontsize',20,'Interpreter','latex')

refer1 = sprintf('Reference %.1f',(convLm1));
refer2 = sprintf('Reference %.1f',(convL2));
legend('L$^{m+1}(\Omega)$-error','L$^{2}(\Omega)$-error',...
    refer1, refer2,'fontsize',16,'Interpreter','latex',...
    'Location','southeast')
refer = sprintf('m = %i',parm);
title(refer,'fontsize',22,'Interpreter','latex')

nnn = sprintf('HDGO%i_m%i_Lmp1',Orden,parm);
export_fig(nnn,'-pdf')

%% CONVERGENCE H1
figure
hold on
plot(Tabla(:,1),Tabla(:,5),'-ob','linewidth',2,'MarkerSize',12)
plot(Tabla(:,1),reference3,'--k','linewidth',2)
set(gca, 'YScale', 'log','XScale', 'log')
set(gca,'fontsize',14)
set(gcf, 'Color', 'w')
box on
axis tight

ax = gca;
ax.YAxis.Exponent = 2;
% NumTicks = 4;
% L = get(gca,'YLim');
% set(gca,'YTick',linspace(ceil(L(1)*10)/10,ceil(L(2)*10)/10,NumTicks))
xtickformat('%.2f')
% ytickformat('%g')
xtickangle(45)

grid on
xlabel('h','fontsize',20,'Interpreter','latex')
ylabel('Error','fontsize',20,'Interpreter','latex')

refer = sprintf('Reference %.1f',(convLH1));
legend('H$^1(\Omega)$-error',refer,'fontsize',16,'Interpreter',...
    'latex','Location','NorthWest')
refer = sprintf('m = %i',parm);
title(refer,'fontsize',22,'Interpreter','latex')

nnn = sprintf('HDGO%i_m%i_H1',Orden,parm);
export_fig(nnn,'-pdf')



