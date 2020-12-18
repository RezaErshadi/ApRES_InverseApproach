clear
close all
clc
[DataPath,RadarData,InvPath,ps] = FUNC_ApRES_PathFix;
tm = string(datetime(now,'ConvertFrom','datenum'));
%%
a = 100;
dB = -a:a;
t = 0:180;
Sr = 1 ./ (tand(t./2)).^2;
Sr_dB = 20.*log10(Sr);

fj = -1.48 .* Sr_dB +90;

pltdim = [0.5000    0.0256    0.3    0.6];
CLASS_FixedPlot.SetFigureSize(pltdim(1),pltdim(2),pltdim(3),pltdim(4));    fg = gcf;   fg.Color = "white";
plot(t,Sr_dB,'b','LineWidth',3)
hold on
plot(fj,Sr_dB,'-r','LineWidth',3)
ylabel('Reflecrion ratio [dB]','Fontsize',20);
xlabel('Nodes angular distance [deg]','Fontsize',20);
fg.InvertHardcopy = 'off';
% grid on
legend('This study','Fujita et al. 2006','location','northeast')
xlim([0 180])
xticks([0:15:180])
set(gca,'FontSize',20)
ylim([-80 80])
yticks([-80:20:80])
set(gca,'FontSize',20)

fg.InvertHardcopy = 'off';
set(fg,'Units','inches');
screenposition = get(fg,'Position');
set(fg,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
print(fg,'ReflecationRatio2AngularDistance','-dpdf','-painters','-r300')
set(gcf,'Units','normalized');
%%
% a = 100;
% dB = -a:a;
% Sr = 10.^(dB./20);
% theta = atan2d(sqrt(Sr),1);
% 
% dt = 2*atan2d(sqrt(Sr),1);
% 
% dBm=(dB(1:end-1)+dB(2:end))/2;
% dtheta = abs(diff(theta));
% 
% fj = 1.4.*dB +45;
% 
% figure,
% fg = gcf;
% fg.Color  = 'white';
% plot(dB,theta,'-b')
% hold on
% plot(dB,fj,'-r')
% xlabel('Scattering anisotropy (dB)','Fontsize',14);
% ylabel('Minimum lateral distance between two nodes (deg)','Fontsize',14);
% fg.InvertHardcopy = 'off';
% grid on
% legend('This study','Fujita et al. 2006','location','northwest')
% ylim([0 90])
% set(gca,'FontSize',14)





















% set(gca,'Color','white')


% % fg.Color = 'black';
% % yyaxis left
% plot(dB,theta,'-b')
% hold on
% plot(dB,fj,'-r')
% % plot(dB,dt,'-k')
% xlabel('Scattering anisotropy (dB)','Color','y','Fontsize',12);
% set(gca,'YColor','y');
% ylabel('\theta (deg)','Color','y','Fontsize',12);
% 
% % yyaxis right
% % plot(dB,dt,'-k')
% % plot(dBm,dtheta)
% % ylabel('\Delta\theta (deg)','Color','r','Fontsize',12);
% % set(gca,'YColor','r');
% fg.InvertHardcopy = 'off';
% set(gca,'XColor','y');
% grid on
% legend('\theta','Fujita et al. 2006','Node Distance')