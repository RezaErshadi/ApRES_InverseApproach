clear;
close all;
clc;
[DataPath,RadarData,InvPath,ps] = FUNC_ApRES_PathFix;
tm = string(datetime(now,'ConvertFrom','datenum'));
%% use forward model
%dZ = 0.1449;
dZ = 1;
% hb = [500,1000,4000];
% Ax = [30;60;90];
% dE = [0.1;0.4;0.3];
% r = [1;1;1];
% hb = [1000];
% Ax = [135]';
% dE = [0.2];
% r = [-1];
% paper synthtic model
hb = [500,1000,1500,2000,2500,3000,4000];
Ax =[45;45;45;45;135;135;120];
dE = [0.025;0.2;0.2;0.2;0.2;0.45;0.2];
r = [0 ; 0 ; 10 ; -10 ; -10 ; -20 ; 0] ;
s2nr =0;
OP = CLASS_FM.BeginForwardModel(hb,dE,r,Ax,dZ,s2nr); 
%%
f = OP.f;
ao = OP.ao;
Z = OP.Z;
dZ = OP.dZ;
Zmx = OP.Zmx;
dtaParams = OP.Dta;
PAHH = dtaParams{5};
PAHV = dtaParams{7};
PCP = dtaParams{14};
Psi = dtaParams{18};
HnodeMin = sort(PAHH(:)); HnodeMin = nanmean(HnodeMin(1:round(0.01*length(HnodeMin),0)));
PAHHnode = PAHH<=HnodeMin;
AxOut = CLASS_Initials.PossiblePrincipalAxis(PAHV ,ao);
% PCA_HV = CLASS_S2P.PCA_denoising(PAHV,1);
% [psbl_Ax2,iax2] = CLASS_Initials.PossiblePrincipalAxis(PCA_HV ,ao);
%%
pltdim = [0.1,0.1,0.8,0.7];
w2plt = [5,14,7,18];
[fg,ax,cb] = CLASS_FixedPlot.AdvancePlot([],dtaParams,ao,Z,w2plt,[],pltdim);

% cla(ax{4})
% imagesc(ax{4},ao,Z,Psi)
caxis(ax{4},[-0.5 0.5])
cm = getPyPlot_cMap('seismic',100);
colormap(ax{4},cm);

cb = colorbar(ax{1});
cb.Label.String = 'P_{HH} [dB]';
cb = colorbar(ax{2});
cb.Label.String = '\phi_{HHVV} [rad]';
cb = colorbar(ax{3});
cb.Label.String = 'P_{HV} [dB]';
cb = colorbar(ax{4});
cb.Label.String = '\Psi [-]';
yticklabels(ax{2},'')
ylabel(ax{2},[])
yticklabels(ax{3},'')
ylabel(ax{3},[])
yticklabels(ax{4},'')
ylabel(ax{4},[])
set(ax{1},'FontSize',18)
set(ax{2},'FontSize',18)
set(ax{3},'FontSize',18)
set(ax{4},'FontSize',18)
title(ax{1},['\fontsize{14}(a) HH Power anomaly'])
title(ax{2},['\fontsize{14}(b) HHVV Coherence phase'])
title(ax{3},['\fontsize{14}(c) HV Power anomaly'])
title(ax{4},['\fontsize{14}(d) Scaled Phase Derivative'])



caxis(ax{1},[HnodeMin -HnodeMin])
caxis(ax{3},[HnodeMin -HnodeMin])
for i = 1:length(hb)
    for j = 1:length(ax)
        hold(ax{j},'on')
        plot(ax{j},[0 180],[hb(i) hb(i)],'--k')
    end
end

% plot(ax{2},AxOut.AxMin(:,1),Z,'.c','MarkerFaceColor','c')
% plot(ax{2},AxOut.AxMin(:,2),Z,'.c','MarkerFaceColor','c')
% plot(ax{3},AxOut.AxMin(:,1),Z,'.c','MarkerFaceColor','c')
% plot(ax{3},AxOut.AxMin(:,2),Z,'.c','MarkerFaceColor','c')
%%
PAHV(isnan(PAHV)) = -inf;
[~,isrt] = sort(PAHV,2);
j1 = 1;
for i = 1:length(hb)
    [~,j2] = min(abs(hb(i)-Z));
    z2 = Z(j1:j2);
    for j = j1:j2
        [~,iao(j,1)] = min(abs(Psi(j,:)-dE(i)));
        dEbest(j,1) = Psi(j,iao(j,1));
    end
    j1=j2+1;
end

ZZ = Z(1):50:Z(end);
j1 = 1;
for i = 1:length(ZZ)
    [~,j2] = min(abs(ZZ(i)-Z));
    psblv(i,1) = isrt(j2,1);
    psblv(i,2) = isrt(j2,2);
    HAbst(i,1) = iao(j2,1);
    j1=j2+1;
end

%%
plot(ax{3},psblv(:,1),ZZ,'sy','MarkerEdgeColor','k','MarkerFaceColor','y','MarkerSize',10)
plot(ax{3},psblv(:,2),ZZ,'sy','MarkerEdgeColor','k','MarkerFaceColor','y','MarkerSize',10)
plot(ax{4},psblv(:,1),ZZ,'sy','MarkerEdgeColor','k','MarkerFaceColor','y','MarkerSize',10)
plot(ax{4},psblv(:,2),ZZ,'sy','MarkerEdgeColor','k','MarkerFaceColor','y','MarkerSize',10)
plot(ax{4},HAbst,ZZ,'ok','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',7)
%%
% dZ2 = 100;
% z2 = 0:dZ2:Z(end);
% i1 = 1;
% for i = 1:length(z2)-1
%     [~,i2] = min(abs(z2(i+1)-Z));
%     sumHV = sum(PAHV(i1:i2,:));
%     [~,imin] = sort(sumHV);
%     nd(i,1) = ao(imin(1));
%     nd(i,2) = ao(imin(2));
%     i1 = i2+1;
% end
% iax = sort(nd,2);
% z2c = z2(1:end-1)+dZ2/2;
% figure,
% plot(nd,z2,'*')
% plot(ax{2},nd,z2c,'dk','MarkerSize',10,'MarkerFaceColor','g')
% plot(ax{4},nd,z2c,'dk','MarkerSize',10,'MarkerFaceColor','g')

fg.InvertHardcopy = 'off';
% print(fg,"fig8.png",'-dpng','-r300')
set(fg,'Units','inches');
screenposition = get(fg,'Position');
set(fg,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
print(fg,'Synthetic7L.pdf','-dpdf','-painters','-r300')
set(gcf,'Units','normalized');
%%
% [otp] = FUNC_MYpca(PAHV,Z,[HnodeMin -HnodeMin],1,0,1,1,1);
%%
% for i = 1:length(Z)
%     psbl_dE1(i,1) = abs(dEhor(i,iax1(i,1)));
%     psbl_dE2(i,1) = abs(dEhor(i,iax2(i,1)));
% end
% psbl_dE1(psbl_dE1>0.5) = nan;
% dE01 = fillmissing(psbl_dE1,'nearest');
% psbl_dE2(psbl_dE2>0.5) = nan;
% dE02 = fillmissing(psbl_dE2,'nearest');
% figure,
% plot(dE02,Z,'sg')
% hold on
% plot(dE01,Z,'.k')
% plot(dEbest,Z,'+b')
% set(gca,'YDir','reverse')
% grid on
% xlim([0 0.5])
% %%
% figure,
% plot(ao,dEhor(1:501,:))
% ylim([0 0.5])
% figure,
% plot(ao,dEhor(502:1001,:))
% ylim([0 0.5])
% figure,
% plot(ao,dEhor(1002:1501,:))
% ylim([0 0.5])
% %%
% dEstd1 = nanstd(dEhor(1:501,:));
% dEstd2 = nanstd(dEhor(502:1001,:));
% dEstd3 = nanstd(dEhor(1002:1501,:));
% figure,
% plot(ao,dEstd1)
% figure,
% plot(ao,dEstd2)
% figure,
% plot(ao,dEstd3)















