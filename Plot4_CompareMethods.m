% This script make the plots of the convergence of HDG vs Others
close all

spr = 'xo*+p><s';
listm = [2 1 4 5 6 7 8];

%% CONVERGENCE LM1
figure
for ij=1:7
    ii = listm(ij);
    Tabla = eval(sprintf('M%i',ii));
    Tabla = Tabla(1:4,:);

    OrdenErrorLm1= abs(log(Tabla(1:(end-1),2)./Tabla(2:(end),2)))'./abs(log(Tabla(1:(end-1),1)./Tabla(2:(end),1)))';

    if ii == 2 
    hold on
    plot(Tabla(:,1),Tabla(:,2),['-',spr(ii),'r'],'linewidth',2,'MarkerSize',12)
    else
    hold on
    plot(Tabla(:,1),Tabla(:,2),['--',spr(ii)],'Color',[0 0 0]+0.04*ii,'linewidth',1.2,'MarkerSize',8)
    end
end

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

legend('HDG','CFVEM','HMM','LEPNC','MLP','VAG1','VAG2',...
    'fontsize',14,'Interpreter','latex',...
    'Location','northeastoutside')
title('L$^{m+1}(\Omega)$-error','fontsize',20,'Interpreter','latex')
 
nnn = sprintf('ComparisonLm1');
export_fig(nnn,'-pdf')

%% CONVERGENCE H1
figure
for ij=1:7
    ii = listm(ij);
    Tabla = eval(sprintf('M%i',ii));
    Tabla = Tabla(1:4,:);

    OrdenErrorH1 = abs(log(Tabla(1:(end-1),3)./Tabla(2:(end),3)))'./abs(log(Tabla(1:(end-1),1)./Tabla(2:(end),1)))';
    
    if ii == 2 
    hold on
    plot(Tabla(:,1),Tabla(:,3),['-',spr(ii),'b'],'linewidth',2,'MarkerSize',12)
    else
    hold on
    plot(Tabla(:,1),Tabla(:,3),['--',spr(ii)],'Color',[0 0 0]+0.04*ii,'linewidth',1.2,'MarkerSize',8)
    end
end

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

legend('HDG','CFVEM','HMM','LEPNC','MLP','VAG1','VAG2',...
    'fontsize',14,'Interpreter','latex',...
    'Location','northeastoutside')
title('H$^{1}(\Omega)$-error','fontsize',20,'Interpreter','latex')

nnn = sprintf('ComparisonH1');
export_fig(nnn,'-pdf')




