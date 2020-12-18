clear;
close all;
clc;
[DataPath,RadarData,InvPath,ps] = FUNC_ApRES_PathFix;
tm = string(datetime(now,'ConvertFrom','datenum'));
%% use forward model
%dZ = 0.1449;
dZ = 1;
% hb = [500,1000];
% Ax = [30;60];
% dE = [0.1;0.4];
% r = [1;1];
% hb = [1000];
% Ax = [150]';
% dE = [0.1];
% r = [-10];
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
% [psbl_Ax1,iax1] = CLASS_Initials.PossiblePrincipalAxis(PAHV ,ao);
% PCA_HV = CLASS_S2P.PCA_denoising(PAHV,1);
% [psbl_Ax2,iax2] = CLASS_Initials.PossiblePrincipalAxis(PCA_HV ,ao);
%%
pltdim = [0.1,0.1,0.5,0.4];
w2plt = [5,14,7,18];
[fg,ax,cb] = CLASS_FixedPlot.AdvancePlot([],dtaParams,ao,Z,w2plt,[],pltdim);
% % for i = 2:length(ax)
% % yticklabels(ax{i},'')
% % ylabel(ax{i},[])
% % end

caxis(ax{1},[HnodeMin -HnodeMin])
caxis(ax{3},[HnodeMin -HnodeMin])

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

set(ax{1},'FontSize',14)
set(ax{2},'FontSize',14)
set(ax{3},'FontSize',14)
set(ax{4},'FontSize',14)

title(ax{1},['\fontsize{14}HH Power anomaly'])
title(ax{3},['\fontsize{14}HV Power anomaly'])
title(ax{2},['\fontsize{14}HHVV Coherence phase'])
title(ax{4},['\fontsize{14}Scaled Phase Derivative'])

for i = 1:length(hb)
hold(ax{1},'on')
plot(ax{1},xlim(ax{1}),[hb(i) hb(i)],'--k','LineWidth',2)
hold(ax{2},'on')
plot(ax{2},xlim(ax{2}),[hb(i) hb(i)],'--k','LineWidth',2)
hold(ax{3},'on')
plot(ax{3},xlim(ax{3}),[hb(i) hb(i)],'--k','LineWidth',2)
hold(ax{4},'on')
plot(ax{4},xlim(ax{4}),[hb(i) hb(i)],'--k','LineWidth',2)
end

% %%
% fg.InvertHardcopy = 'off';
% print(fg,"Synthbased(r=10,v1=135).png",'-dpng','-r300')






