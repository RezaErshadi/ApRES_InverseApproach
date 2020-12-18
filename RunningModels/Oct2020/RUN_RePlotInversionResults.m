clear
close all
clc
[DataPath,RadarData,InvPath,ps] = FUNC_ApRES_PathFix;
tm = string(datetime(now,'ConvertFrom','datenum'));
%%
% temp2(1).name = 'PpRES_Cores_PpRES_LD01_12-Nov-2020 17:30:30_InversionResults.mat';
% temp3 = strcat(InvPath,ps,'Corrected_HHVV',ps,'EDML_FinalInversion');
% FD = [80 20]; % flow direction

temp2(1).name = '11_Concordia_BAS_EPICA_12-Nov-2020 17:45:02_InversionResults.mat';
temp3 = strcat(InvPath,ps,'Corrected_HHVV',ps,'DC_FinalInversion');
FD = [35 20]; % flow direction

% temp2 = dir(strcat(temp3,ps,'Data',ps,'*.mat'));
%%
for i = 1:size(temp2,1)
FileName = temp2(i).name;
FilePath = strcat(temp3,ps,'Data',ps,FileName);
%%
fp = FilePath;
% fp = strcat(InvPath,ps,FilePath);
load(fp);
%% Read data
%---------------------------------------------------------------------------------------------
% MetaData of the ApRES data
MetaData = InversionOutput.MetaData;
%---------------------------------------------------------------------------------------------
% All observed parameters
OP_ObserevedData = InversionOutput.OP_ObserevedData;
obs_f = InversionOutput.OP_ObserevedData.f; % frequency
obs_ao = InversionOutput.OP_ObserevedData.ao; % azimuth
obs_Z = InversionOutput.OP_ObserevedData.Z; % depth
obs_dz = InversionOutput.OP_ObserevedData.dZ; % depth spacing
obs_Zmx = InversionOutput.OP_ObserevedData.Zmx; % max depth
% Signal parameters including: 
% HH,VH,HV,VV,PAHH,PAVH,PAHV,PAVV,PDHH,PDVH,PDHV,PDVV,absC,argC,reC,imC,gradC,Psi
obs_Dta = InversionOutput.OP_ObserevedData.Dta; 
obsPAHH = obs_Dta{5};
obsPAHV = obs_Dta{7};
obsCP = obs_Dta{14};
obsPsi = obs_Dta{18};
%---------------------------------------------------------------------------------------------
% All modeled parameters
OP_OptimizedModel = InversionOutput.OP_OptimizedModel;
mdl_f = InversionOutput.OP_OptimizedModel.f; % frequency
mdl_ao = InversionOutput.OP_OptimizedModel.ao; % azimuth
mdl_Z = InversionOutput.OP_OptimizedModel.Z; % depth
mdl_dz = InversionOutput.OP_OptimizedModel.dZ; % depth spacing
mdl_Zmx = InversionOutput.OP_OptimizedModel.Zmx; % max depth
mdl_p1 = InversionOutput.OP_OptimizedModel.P1; % p1 factor
mdl_p2 = InversionOutput.OP_OptimizedModel.P2; % p2 factor
mdl_mp = InversionOutput.OP_OptimizedModel.mp; % model parameters
% Signal parameters including: 
% HH,VH,HV,VV,PAHH,PAVH,PAHV,PAVV,PDHH,PDVH,PDHV,PDVV,absC,argC,reC,imC,gradC,Psi
mdl_Dta = InversionOutput.OP_OptimizedModel.Dta; % Signal parameters
mdlPAHH = mdl_Dta{5};
mdlPAHV = mdl_Dta{7};
mdlCP = mdl_Dta{14};
mdlPsi = mdl_Dta{18};
%---------------------------------------------------------------------------------------------
coreEV = InversionOutput.coreEV; % Ice core EigenValues
%---------------------------------------------------------------------------------------------
Zmdl = InversionOutput.Zmdl; % model depth vector
AxOut = InversionOutput.AxOut; % min of HV power anomaly

HA0 = InversionOutput.HA0; % Initial Horizontal Anisotropy (lambda2-lambda1) at min of HV power anomaly
r0 = InversionOutput.r0; % Initial reflection ratio (r)
v1_0 = InversionOutput.v1_0; % Secondary guess for first Eigenvector orientation (v1)

r00 = InversionOutput.r00; %  optimized reflection ratio from inverse approach
HA00 = InversionOutput.HA00; %  optimized horizontal anisotropy from inverse approach

HA = InversionOutput.HA; % final horizontal anisotropy
r = InversionOutput.r; % final hreflection ratio
v1 = InversionOutput.v1; % optimized first Eigenvector orientation from inverse approach (v1)
v2 = InversionOutput.v2; % final estimation for second Eigenvector orientation (v2)

EigVal = InversionOutput.EigVal; % estimated Eigenvalues
%% PLOT
[fg,ax] = CLASS_InvPlot.InversionPlot(obsPAHH,obsPAHV,obsCP,obsPsi,obs_Z,obs_ao,[],1);

ax = CLASS_InvPlot.UpdateInversionPlot(ax,mdlPAHH,mdlPAHV,mdlCP,mdlPsi,mdl_Z,[],mdl_ao,[],[],1);

p1 = plot(ax{3},AxOut.AxMin(:,1),obs_Z,'.','color',[0.309, 1, 0.039],'MarkerFaceColor',[0.309, 1, 0.039]);
p2 = plot(ax{3},AxOut.AxMin(:,2),obs_Z,'.','color',[0.309, 1, 0.039],'MarkerFaceColor',[0.309, 1, 0.039]);
p3 = plot(ax{4},AxOut.AxMin(:,1),obs_Z,'.','color',[0.309, 1, 0.039],'MarkerFaceColor',[0.309, 1, 0.039]);
p4 = plot(ax{4},AxOut.AxMin(:,2),obs_Z,'.','color',[0.309, 1, 0.039],'MarkerFaceColor',[0.309, 1, 0.039]);

p5 = plot(ax{9},coreEV(:,2),coreEV(:,1),'xk');
p6 = plot(ax{9},coreEV(:,3),coreEV(:,1),'dk');
p7 = plot(ax{9},coreEV(:,4),coreEV(:,1),'*k');

p8 = plot(ax{9},EigVal(:,1),Zmdl,'.-r','MarkerFaceColor','r','LineWidth',2,'MarkerSize',15);
p9 = plot(ax{9},EigVal(:,2),Zmdl,'.-b','MarkerFaceColor','b','LineWidth',2,'MarkerSize',15);
p10 = plot(ax{9},EigVal(:,3),Zmdl,'.-g','MarkerFaceColor','g','LineWidth',2,'MarkerSize',15);

p11 = plot(ax{10},coreEV(:,5),coreEV(:,1),'xk');
p12 = plot(ax{10},HA,Zmdl,'.-r','LineWidth',2,'MarkerSize',15);

p13 = rectangle(ax{11},'Position',[FD(1) obs_Z(1) FD(2) obs_Z(end)],'FaceColor',[0.984, 0.6, 0.054 0.75],'EdgeColor','k');

p13legend = line(NaN,NaN,'LineWidth',5,'Color',[0.984, 0.6, 0.054 1]);

p14 = plot(ax{11},v1,Zmdl,'.-b','LineWidth',2,'MarkerSize',15);
p15 = plot(ax{11},v2,Zmdl,'.-r','LineWidth',2,'MarkerSize',15);

p16 = plot(ax{12},r,Zmdl,'.-r','LineWidth',2,'MarkerSize',15);

lgnbck = [0.9 0.9 0.9];
legend(ax{9},[p5 p6 p7 p8 p9 p10],...
    {'Meas. \lambda1','Meas. \lambda2','Meas. \lambda3','Est. \lambda1','Est. \lambda2','Est. \lambda3'},...
    'FontSize',11,'FontWeight','bold',...
    'Position',[0.24571250368958 0.259042995071674 0.0616432424968738 0.0666436485840816],'Color',lgnbck)
legend(ax{10},[p11 p12],{'Measured','Estimated'},'FontSize',11,'FontWeight','bold',...
    'Position',[0.425078331806832 0.298086246235756 0.0652871780746877 0.0279696141297652],'Color',lgnbck)
legend(ax{11},[p14 p15 p13legend],{'Est. v1','Est. v2' 'Flow dir.'},'FontSize',11,'FontWeight','bold',...
    'Position',[0.672292278281226 0.291352820712628 0.0606021180460697 0.0345303878279524],'Color',lgnbck)

fg.InvertHardcopy = 'off';
%%
fp = fp(1:end-4);
fp = replace(fp,'Data','Figures');
set(fg,'Units','inches');
screenposition = get(fg,'Position');
set(fg,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
print(fg,fp,'-dpdf','-painters','-r300')
set(gcf,'Units','normalized');
% close all
end



