% This script make the plots of the Negative mass 
spr = 'bmg';

figure
for aa= 1:3
    rrr = sprintf('Tab_O1_m%i.mat',aa+1);
    load(rrr)
    
    hold on
    plot(Tabla(:,1),Tabla(:,8),['-*',spr(aa)],'linewidth',2,'MarkerSize',12)
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
ylabel('NMass','fontsize',20,'Interpreter','latex')

legend('$m=2$','$m=3$','$m=4$','fontsize',18,'Interpreter',...
    'latex','Location','NorthWest')

nnn = sprintf('HDGO%i_Neg',Orden);
export_fig(nnn,'-pdf')
