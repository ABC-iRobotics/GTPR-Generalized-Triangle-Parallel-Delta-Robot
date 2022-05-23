
% Example1_optimized;
% Example2_optimized;

load("case1.mat");
load("case2.mat");



% figure(6);
% xticks(xAxis_main);
% yticks(yAxis_main);
% zticks(zAxis_main);
% set(gca, 'XTickLabel', string(xAxis_main));
% set(gca, 'YTickLabel', string(yAxis_main));
% set(gca, 'ZTickLabel', string(zAxis_main));
% xtickformat('%.1f');
% ytickformat('%.1f');
% ztickformat('%.1f');
% set(gcf,'Color','w');
% saveas(gcf, 'case1_fig6.eps', 'epsc');
% saveas(gcf, 'case1_fig6.png');




% 
% figure(7);
% xticks(xAxis_main);
% yticks(yAxis_main);
% zticks(zAxis_main);
% set(gca, 'XTickLabel', string(xAxis_main));
% set(gca, 'YTickLabel', string(yAxis_main));
% set(gca, 'ZTickLabel', string(zAxis_main));
% xtickformat('%.1f');
% ytickformat('%.1f');
% ztickformat('%.1f');
% set(gcf,'Color','w');
% saveas(gcf, 'case2_fig1.eps', 'epsc');
% saveas(gcf, 'case2_fig1.png');
% 
% figure(8);
% xticks(xAxis_main);
% yticks(yAxis_main);
% zticks(zAxis_main);
% set(gca, 'XTickLabel', string(xAxis_main));
% set(gca, 'YTickLabel', string(yAxis_main));
% set(gca, 'ZTickLabel', string(zAxis_main));
% xtickformat('%.1f');
% ytickformat('%.1f');
% ztickformat('%.1f');
% set(gcf,'Color','w');
% saveas(gcf, 'case2_fig2.eps', 'epsc');
% saveas(gcf, 'case2_fig2.png');
% 
% 
% figure(9);
% xticks(xAxis_main);
% yticks(yAxis_main);
% zticks(zAxis_main);
% set(gca, 'XTickLabel', string(xAxis_main));
% set(gca, 'YTickLabel', string(yAxis_main));
% set(gca, 'ZTickLabel', string(zAxis_main));
% xtickformat('%.1f');
% ytickformat('%.1f');
% ztickformat('%.1f');
% set(gcf,'Color','w');
% saveas(gcf, 'case2_fig3.eps', 'epsc');
% saveas(gcf, 'case2_fig3.png');
% 
% 
% figure(10);
% xticks(xAxis_main);
% yticks(yAxis_main);
% zticks(zAxis_main);
% set(gca, 'XTickLabel', string(xAxis_main));
% set(gca, 'YTickLabel', string(yAxis_main));
% set(gca, 'ZTickLabel', string(zAxis_main));
% xtickformat('%.1f');
% ytickformat('%.1f');
% ztickformat('%.1f');
% set(gcf,'Color','w');
% saveas(gcf, 'case2_fig4.eps', 'epsc');
% saveas(gcf, 'case2_fig4.png');
% 
% 
% figure(11);
% xticks(xAxis_main);
% yticks(yAxis_main);
% zticks(zAxis_main);
% set(gca, 'XTickLabel', string(xAxis_main));
% set(gca, 'YTickLabel', string(yAxis_main));
% set(gca, 'ZTickLabel', string(zAxis_main));
% xtickformat('%.1f');
% ytickformat('%.1f');
% ztickformat('%.1f');
% set(gcf,'Color','w');
% saveas(gcf, 'case2_fig5.eps', 'epsc');
% saveas(gcf, 'case2_fig5.png');
% 
% 
% figure(12);
% xticks(xAxis_main);
% yticks(yAxis_main);
% zticks(zAxis_main);
% set(gca, 'XTickLabel', string(xAxis_main));
% set(gca, 'YTickLabel', string(yAxis_main));
% set(gca, 'ZTickLabel', string(zAxis_main));
% xtickformat('%.1f');
% ytickformat('%.1f');
% ztickformat('%.1f');
% set(gcf,'Color','w');
% saveas(gcf, 'case2_fig6.eps', 'epsc');
% saveas(gcf, 'case2_fig6.png');

figure(6);
xAxis1Plot = xticks;
yAxis1Plot = yticks;
zAxis1Plot = zticks;


figure(12);
xAxis2Plot = round(xticks,1);
yAxis2Plot = round(yticks,1);
zAxis2Plot = round(zticks,1);

% figure(1);
% xlim([min([xAxis1Plot, xAxis2Plot]) max([xAxis1Plot, xAxis2Plot])])
% ylim([min([yAxis1Plot, yAxis2Plot]) max([yAxis1Plot, yAxis2Plot])])
% zlim([min([zAxis1Plot, zAxis2Plot]) max([zAxis1Plot, zAxis2Plot])])

% xticks('auto');
% yticks('auto');
% zticks('auto');

% xtickformat('%.1f');
% ytickformat('%.1f');
% ztickformat('%.1f');
% set(gcf,'Color','w');
% saveas(gcf, 'case1_fig1.eps', 'epsc');
% saveas(gcf, 'case1_fig1.png');

% limsx = [min([xAxis1Plot, xAxis2Plot])-0.1 max([xAxis1Plot, xAxis2Plot])+0.1];
% limsy = [min([yAxis1Plot, yAxis2Plot])-0.1 max([yAxis1Plot, yAxis2Plot])+0.1];
% limsz = [min([zAxis1Plot, zAxis2Plot])-0.1 max([zAxis1Plot, zAxis2Plot])+0.1];

limsx = [min([xAxis1, xAxis2]) max([xAxis1, xAxis2])];
limsy = [min([yAxis1, yAxis2]) max([yAxis1, yAxis2])];
limsz = [min([zAxis1, zAxis2]) max([zAxis1, zAxis2])];

ticksx = limsx(1):0.1:limsx(2);
ticksy = limsy(1):0.1:limsy(2);
ticksz = limsz(1):0.1:limsz(2);



%%
figure(4);
xlim(limsx)
ylim(limsy)
zlim(limsz)

xticks(ticksx);
yticks(ticksy);
zticks(ticksz);

% xticks('auto');
% yticks('auto');
% zticks('auto');

xtickformat('%.1f');
ytickformat('%.1f');
ztickformat('%.1f');
set(gcf,'Color','w');
saveas(gcf, 'case1_fig4.eps', 'epsc');
saveas(gcf, 'case1_fig4.png');

%%
figure(5);
xlim(limsx)
ylim(limsy)
zlim(limsz)

xticks(ticksx);
yticks(ticksy);
zticks(ticksz);
xtickformat('%.1f');
ytickformat('%.1f');
ztickformat('%.1f');
set(gcf,'Color','w');
saveas(gcf, 'case1_fig5.eps', 'epsc');
saveas(gcf, 'case2_fig5.png');


%%
figure(6);
xlim(limsx)
ylim(limsy)
zlim(limsz)

xticks(ticksx);
yticks(ticksy);
zticks(ticksz);

xtickformat('%.1f');
ytickformat('%.1f');
ztickformat('%.1f');
set(gcf,'Color','w');
saveas(gcf, 'case1_fig6.eps', 'epsc');
saveas(gcf, 'case1_fig6.png');

%%
figure(10);
xlim(limsx)
ylim(limsy)
zlim(limsz)

xticks(ticksx);
yticks(ticksy);
zticks(ticksz);

xtickformat('%.1f');
ytickformat('%.1f');
ztickformat('%.1f');
set(gcf,'Color','w');
saveas(gcf, 'case2_fig10.eps', 'epsc');
saveas(gcf, 'case2_fig10.png');

%%
figure(11);
xlim(limsx)
ylim(limsy)
zlim(limsz)

xticks(ticksx);
yticks(ticksy);
zticks(ticksz);

xtickformat('%.1f');
ytickformat('%.1f');
ztickformat('%.1f');
set(gcf,'Color','w');
saveas(gcf, 'case2_fig11.eps', 'epsc');
saveas(gcf, 'case2_fig11.png');

%%
figure(12);
xlim(limsx)
ylim(limsy)
zlim(limsz)

xticks(ticksx);
yticks(ticksy);
zticks(ticksz);

xtickformat('%.1f');
ytickformat('%.1f');
ztickformat('%.1f');
set(gcf,'Color','w');
saveas(gcf, 'case2_fig12.eps', 'epsc');
saveas(gcf, 'case2_fig12.png');


%%

figure(1);
xtickformat('%.1f');
ytickformat('%.1f');
ztickformat('%.1f');
set(gcf,'Color','w');
xlim(limsx)
ylim(limsy)
zlim(limsz)

xticks(ticksx);
yticks(ticksy);
zticks(ticksz);
saveas(gcf, 'case1_fig1.eps', 'epsc');
saveas(gcf, 'case1_fig1.png');

figure(7);
xlim(limsx)
ylim(limsy)
zlim(limsz)

xticks(ticksx);
yticks(ticksy);
zticks(ticksz);
xtickformat('%.1f');
ytickformat('%.1f');
ztickformat('%.1f');
set(gcf,'Color','w');
saveas(gcf, 'case2_fig7.eps', 'epsc');
saveas(gcf, 'case2_fig7.png');






