clear
close all
clc
[DataPath,RadarData,InvPath,ps] = FUNC_ApRES_PathFix;
tm = string(datetime(now,'ConvertFrom','datenum'));
%%
load('E1E2E3_Schmidt.mat')
% E2E1 = (0:0.001:0.5);
% E3E2 = (0:0.001:1);
% [X,Y] = meshgrid(E2E1,E3E2);
% for i = 1:length(E2E1)
%     tic
%     parfor j = 1:length(E3E2)
%         x = E2E1(i);
%         y = E3E2(j);
%         e3 = (1+(2*(y))+x)/3;
%         e2 = e3-y;
%         e1 = e2-x;
%         check1 = (e1 < 0 || e2 < 0 || e3 < 0);
%         check2 = (e1 > 1 || e2 > 1 || e3 > 1);
%         check3 = (e3<e1);
%         check4 = (e3<e2);
%         check5 = (e2<e1);
%         sumcheck = sum([check1,check2,check3,check4,check5]);
%         if sumcheck ~= 0
%             e1 = nan;
%             e2 = nan;
%             e3 = nan;
%         end
%         E1(j,i) = e1;
%         E2(j,i) = e2;
%         E3(j,i) = e3;
%     end
%     toc
% end
%%
cm = getPyPlot_cMap('seismic',21);
figure,
fg = gcf;
fg.Color = "white";
set(fg,'position',[1922         822        1918         694])
ax1 = subplot(1,3,1);
s = pcolor(ax1,X,Y,E1);
set(s, 'EdgeColor', 'none');
xlabel('\lambda2-\lambda1')
ylabel('\lambda3-\lambda2')
cb1 = colorbar;
cb1.Ticks = [min(E1(:)),max(E1(:))];
% cb1.TickLabels = {string(round(cb1.Ticks(1),2))+" (min E1)",string(round(cb1.Ticks(2),2)),string(round(cb1.Ticks(3),2))+" (max E1)"};
cb1.TickLabels = {string(round(cb1.Ticks(1),2))+" (min \lambda1)",string(round(cb1.Ticks(2),2))+" (max \lambda1)"};
colormap(cm)
caxis([0 1])
title('(a) \lambda1')
set(ax1,'FontSize',16)
pbaspect([1 2 1])
% hold on
% plot(0,0,'ok','MarkerSize',20,'MarkerFaceColor',[0.5 0.5 0.5]);
% plot(0,1,'ok','MarkerSize',20,'MarkerFaceColor',[0.5 0.5 0.5]);
% plot(0.5,0,'ok','MarkerSize',20,'MarkerFaceColor',[0.5 0.5 0.5]);
xticks(ax1,[0:0.1:0.5]);
yticks(ax1,[0:0.1:1]);
set(ax1,'color',[0.9 0.9 0.9])
grid(ax1,'on')
set(ax1,'layer','top')

ax2 = subplot(1,3,2);
s = pcolor(ax2,X,Y,E2);
set(s, 'EdgeColor', 'none');
xlabel('\lambda2-\lambda1')
% ylabel('E3-E2')
cb2 = colorbar;
cb2.Ticks = [min(E2(:)),max(E2(:))];
cb2.TickLabels = {string(round(cb2.Ticks(1),2))+" (min \lambda2)",string(round(cb2.Ticks(2),2))+" (max \lambda2)"};
colormap(cm)
caxis([0 1])
title('(b) \lambda2')
set(ax2,'FontSize',16)
xticks(ax2,[0:0.1:0.5]);
yticks(ax2,[0:0.2:1]);
yticklabels(ax2,[]);
set(ax2,'color',[0.9 0.9 0.9])
grid(ax2,'on')
pbaspect([1 2 1])
set(ax2,'layer','top')

ax3 = subplot(1,3,3);
s = pcolor(ax3,X,Y,E3);
set(s, 'EdgeColor', 'none');
xlabel('\lambda2-\lambda1')
% ylabel('E3-E2')
cb3 = colorbar;
cb3.Ticks = [0 0.33 0.5 1];
% cb3.Label.String = "Eigenvalue";
cb3.Ticks = [min(E3(:)),max(E3(:))];
cb3.TickLabels = {string(round(cb3.Ticks(1),2))+" (min \lambda3)",string(round(cb3.Ticks(2),2))+" (max \lambda3)"};
colormap(cm)
caxis([0 1])
title('(c) \lambda3')
set(ax3,'FontSize',16)
xticks(ax3,[0:0.1:0.5]);
yticks(ax3,[0:0.2:1]);
yticklabels(ax3,[]);
set(ax3,'color',[0.9 0.9 0.9])
grid(ax3,'on')
pbaspect([1 2 1])
set(ax3,'layer','top')


% annotation(fg,'textbox',[0.137448132780083 0.098853868194842 0.0212655601659751 0.0606045845272198],'Color',[1 1 1],...
%     'String',{'A'},...
%     'FontWeight','bold',...
%     'FontSize',20,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');
% 
% annotation(fg,'textbox',[0.137966804979253 0.872492836676217 0.0212655601659751 0.0606045845272197],'Color',[1 1 1],...
%     'String','B',...
%     'FontWeight','bold',...
%     'FontSize',20,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');
% 
% annotation(fg,'textbox',[0.271784232365145 0.103151862464182 0.021265560165975 0.0606045845272197],'Color',[1 1 1],...
%     'String','C',...
%     'FontWeight','bold',...
%     'FontSize',20,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');


% hold(ax1,'on')
% plot(ax1,0,0,'sr','MarkerSize',15,'MarkerFaceColor','r')
% plot(ax1,0,1,'dr','MarkerSize',15,'MarkerFaceColor','r')
% plot(ax1,0.5,0,'^','MarkerSize',15,'MarkerFaceColor','r')

%%
fg.InvertHardcopy = 'off';

% set(fg,'Units','inches');
% screenposition = get(fg,'Position');
% set(fg,...
%         'PaperPosition',[0 0 screenposition(3:4)],...
%         'PaperSize',[screenposition(3:4)]);
% print(fg,'FabricType','-dpdf','-painters','-r300')
% set(gcf,'Units','normalized');

print(fg,"FabricType.png",'-dpng','-r300')







