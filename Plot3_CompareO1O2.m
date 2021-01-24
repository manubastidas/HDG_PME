% This script make the plots of the comparison order 1 and order 2
close all 

load('Tab_O1_m4.mat')
TablaO1 = Tabla;
load('Tab_O2_m4.mat')
TablaO2 = Tabla;

%% Errors O1
OrdenErrorLm1_O1= abs(log(TablaO1(1:(end-1),3)./TablaO1(2:(end),3)))'./abs(log(TablaO1(1:(end-1),1)./TablaO1(2:(end),1)))';
OrdenErrorH1_O1 = abs(log(TablaO1(1:(end-1),5)./TablaO1(2:(end),5)))'./abs(log(TablaO1(1:(end-1),1)./TablaO1(2:(end),1)))';

convLm1_O1 = ceil(max(OrdenErrorLm1_O1(end))*10)/10;
reference1 = TablaO1(:,1).^(convLm1_O1);

convLH1_O1 = ceil(max(OrdenErrorH1_O1(end))*10)/10;
reference2 = TablaO1(:,1.^(convLH1_O1));

%% Errors O2
OrdenErrorLm1_O2= abs(log(TablaO2(1:(end-1),3)./TablaO2(2:(end),3)))'./abs(log(TablaO2(1:(end-1),1)./TablaO2(2:(end),1)))';
OrdenErrorH1_O2 = abs(log(TablaO2(1:(end-1),5)./TablaO2(2:(end),5)))'./abs(log(TablaO2(1:(end-1),1)./TablaO2(2:(end),1)))';

convLm1_O2  = ceil(max(OrdenErrorLm1_O2(end))*10)/10;
reference11 = TablaO2(:,1).^(convLm1_O2);

convLH1_O2  = ceil(max(OrdenErrorH1_O2(end))*10)/10;
reference22 = TablaO2(:,1).^(convLH1_O2);

%% CONVERGENCE LM+1
figure
hold on
plot(TablaO1(:,1),TablaO1(:,3),'-sb','linewidth',2,'MarkerSize',12)
plot(TablaO2(:,1),TablaO2(:,3),'-or','linewidth',2,'MarkerSize',12)
plot(TablaO1(:,1),reference1,'-.k','linewidth',2)
% plot(TablaO2(:,1),reference11,'--k','linewidth',2)
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

refer1 = sprintf('Ref %.1f',(convLm1_O1));
% refer2 = sprintf('Ref %.1f (r=2)',(convLm1_O2));
legend('L$^{m+1}(\Omega)$-error (r=1)','L$^{m+1}(\Omega)$-error (r=2)',...
    refer1,'fontsize',16,'Interpreter','latex',...
    'Location','southeast')
refer = sprintf('m = %i',parm);
title(refer,'fontsize',22,'Interpreter','latex')

nnn = sprintf('HDGO%i_m%iO1O2_Lmp1',Orden,parm);
export_fig(nnn,'-pdf')

%% CONVERGENCE H1
figure
hold on
plot(TablaO1(:,1),TablaO1(:,5),'-sb','linewidth',2,'MarkerSize',12)
plot(TablaO2(:,1),TablaO2(:,5),'-or','linewidth',2,'MarkerSize',12)
plot(TablaO1(:,1),reference2,'-.k','linewidth',2)
plot(TablaO2(:,1),reference22,'--k','linewidth',2)
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

refer1 = sprintf('Ref %.1f (r=1)',(convLH1_O1));
refer2 = sprintf('Ref %.1f (r=2)',(convLH1_O2));
legend('H$^1(\Omega)$-error (r=1)','H$^1(\Omega)$-error (r=2)',...
    refer1, refer2,'fontsize',16,'Interpreter','latex',...
    'Location','southeast')
refer = sprintf('m = %i',parm);
title(refer,'fontsize',22,'Interpreter','latex')

nnn = sprintf('HDGO%i_m%iO1O2_H1',Orden,parm);
export_fig(nnn,'-pdf')