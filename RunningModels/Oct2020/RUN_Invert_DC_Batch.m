clear
close all
clc
[DataPath,RadarData,InvPath,ps] = FUNC_ApRES_PathFix;
tm = string(datetime(now,'ConvertFrom','datenum'));
%% -------------------------------- Data Path
% Path of ApRES data
mainDataPath = '/mnt/esd01/docs/mershadi/WorkDirectory/ApRES/Data/RadarData_ApRES/Concordia_BAS/';
% mainDataPath = '/home/reza/OneDrive(bwedu)/WorkDirectory/ApRES/Data/RadarData_ApRES/Concordia_BAS/';
ff = ["W18","W12","W06","W4d5","W2d5","W1d5","W1d0","W0d5","E0","EPICA","E0d5","E1d0","E1d5","E02","E03","E4d5","E06","E09","E12","E18"];
ss = ["1","2","3","4","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21"];
%% -------------------------------- Input Parameters
maxRange = 2001;
% FigVis = 1;
% kk = 16;
% iss = find(str2double(ss) == kk);
% DtaDir = strcat(mainDataPath,ff(iss));
% FUNC_Batch_invert(DtaDir,maxRange,FigVis,ps,InvPath,tm);

% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print -dpdf -painters epsFig
% set(gcf,'Units','normalized');

FigVis = 0;
for kk = 1:length(ss)
    DtaDir = strcat(mainDataPath,ff(kk));
    SiteNumber = str2double(ss(kk));
    FUNC_Batch_invert(DtaDir,maxRange,FigVis,ps,InvPath,tm);
end
%%
function [] = FUNC_Batch_invert(DtaDir,maxRange,FigVis,ps,InvPath,tm)
EvDir = 'Data/EigenValues/EPICADomeC_EigenValues_Maurine.txt';
BedRange = [];
ResZmdl = 50; % Model Depth Resolution
dA = 1;
ao = 0:dA:179; 
C_DepthWin = maxRange * 0.05;
C_ConvWin = maxRange * 0.05;
DenoisingFlag.PA = ["1", "MovingAverage" , string(maxRange*0.05) ;
                                      "0", "Conv1D" , string(maxRange*0.01) ;
                                      "2", "Conv2D" , string(maxRange*0.05) ;
                                      "0", "DenoisePCA" , string(1)];
DenoisingFlag.PD = ["1", "MovingAverage" , string(maxRange*0.05) ;
                                      "0", "Conv1D" , string(maxRange*0.01) ;
                                      "0", "Conv2D" , string(maxRange*0.01) ;
                                      "0", "DenoisePCA" , string(1)];
%% -------------------------------- Read the Data
[hh,vh,hv,vv,Z,dZ,maz,f,Bed,MetaData] = FUNC_ReadApRES(DtaDir,maxRange,BedRange);
%-------------------
[~,izmx] = min(abs(Z-maxRange));
Z = Z(1:izmx,:); 
hh = hh(1:izmx,:); 
vh = hv(1:izmx,:); 
hv = hv(1:izmx,:); 
vv = vv(1:izmx,:);
SiteNumber = str2double(MetaData(MetaData(:,1) == "SiteNumber",2));
HH0 = hh(:,maz==0); 
VH0 = vh(:,maz==0);
HV0 = hv(:,maz==0); 
VV0 = vv(:,maz==0); 
GeoRef = str2double(MetaData(MetaData(:,1) == "TR2TN",2));
GeoRef = GeoRef +90;
[HH,VV,HV,VH] = CLASS_S2P.AzimuthSynthesizer(HH0,VH0,HV0,VV0,ao,GeoRef);
ObsDta = CLASS_S2P.Signal2Param(HH,VH,HV,VV,Z,ao,f,C_DepthWin,C_ConvWin,DenoisingFlag,"radar");
RadarDta.f = f;
RadarDta.ao = ao;
RadarDta.Z = Z;
RadarDta.dZ = dZ;
RadarDta.Zmx = max(Z);
RadarDta.Dta = ObsDta;
Zmx = max(Z);
PAHH = ObsDta{5};
PAHV = ObsDta{7};
CP = ObsDta{14};
Psi = ObsDta{18};
FoldData = split(DtaDir,ps);
LabelName = strcat(InvPath,ps,string(SiteNumber),'_',FoldData{end-1},'_',FoldData{end},'_',tm);
%%
cp1 = CP;
aomat = repmat(ao,length(Z),1);
Zmat = repmat(Z,1,length(ao));
[~,cp2] = min(abs(pi-cp1(:)));
NodeInfo(1,:) = aomat(cp2);
NodeInfo(2,:) = Zmat(cp2);
ZNode = mean(NodeInfo(2,:));
cp3 = cp1(1:end-1,:).*cp1(2:end,:);
cp3 = [cp1(1,:) ; cp3];
cp4 = find(cp3<0);
BorderInfo(1,:) = aomat(cp4);
BorderInfo(2,:) = Zmat(cp4);
cp5 = zeros(size(CP));
cp5(cp4) = 1;
cp5 = [Z sum(cp5,2)];
cp5(cp5(:,1)<100 | cp5(:,1)>500,:) = [];
[~,imhb] = max(cp5(:,2));
Ztop = cp5(imhb,1);
%% --------------------------------
EV = CLASS_ApRES_Data.LoadEigenValues(EvDir);
% MHB = funcHBDC(SiteNumber,Zmx);
% MHB = [200;400;600; Zmx];
MHB = [Ztop ; Zmx];
Zmdl = (ResZmdl:ResZmdl:Zmx)'; Zmdl(end) = Zmx;
% --------------------------------
[fg,ax] = CLASS_InvPlot.InversionPlot(PAHH,PAHV,CP,Psi,Z,ao,[],FigVis);
% plot ice core Eigenvalues
plot(ax{9},EV(:,2),EV(:,1),'-sk','MarkerFaceColor','y');
plot(ax{9},EV(:,3),EV(:,1),'-sk','MarkerFaceColor','c');
plot(ax{9},EV(:,4),EV(:,1),'-sk','MarkerFaceColor','m');
% plot ice core lambda2-lambda1
plot(ax{10},EV(:,5),EV(:,1),'-sk','MarkerFaceColor','w');
AxOut = CLASS_Initials.PossiblePrincipalAxis(PAHV,ao);
% plot min of HV power anomaly on: dPHV and Psi
plot(ax{2},BorderInfo(1,:),BorderInfo(2,:),'.k','linewidth',2)
plot(ax{1},linspace(ao(1),ao(end),20) , ZNode ,'*w','linewidth',2)
plot(ax{2},linspace(ao(1),ao(end),20) , ZNode ,'*w','linewidth',2)
plot(ax{3},AxOut.AxMin(:,1),Z,'.c','MarkerFaceColor','c')
plot(ax{3},AxOut.AxMin(:,2),Z,'.c','MarkerFaceColor','c')
plot(ax{4},AxOut.AxMin(:,1),Z,'.c','MarkerFaceColor','c')
plot(ax{4},AxOut.AxMin(:,2),Z,'.c','MarkerFaceColor','c')
%--------------------------------- Add North and surface velocity direction
SV = [35 20];
rectangle(ax{11},'Position',[SV(1) Z(1) SV(2) Z(end)],'FaceColor',[1 1 0 0.4],'EdgeColor','k')
% --------------------------------
HA0 = CLASS_Initials.PossibledE(Z,Zmdl,AxOut,Psi); % Initial Horizontal Anisotropy
% Plot initial HA
plot(ax{10},HA0,Zmdl,'sb','MarkerFaceColor','b'); drawnow 
% --------------------------------
r0 = zeros(size(Zmdl)); % Initial Reflection Ratio
% Plot initial r
plot(ax{12},r0,Zmdl,'sb','MarkerFaceColor','b'); drawnow 
% --------------------------------
[v1_0,ax] = FUNC_SpecificLayers(MHB,HA0,r0,Zmdl,ObsDta,dZ,ax,[14]); % Initial Eigenvector
% Plot initial v1
plot(ax{11},v1_0,Zmdl,'sb','MarkerFaceColor','b'); drawnow
% --------------------------------
ax = FUNC_PlotOptimizedModel(ax,HA0,r0,v1_0,Zmdl,Z,dZ,ao);
% -------------------------------- INITIALIZING INVERRSION
options1 = optimoptions('fmincon','Display','iter-detailed');
options1.UseParallel = true;
options1.OptimalityTolerance = 1e-6;
options1.FunctionTolerance = 1e-6;
options1.ConstraintTolerance = 1e12;
%-------------------------------- INVERSION
resInv = ResZmdl;
Zinv = (resInv:resInv:Z(end))';
%% -------------------------------- INVERT v1 & v2
if FigVis == 1
    cla(ax{11})
    rectangle(ax{11},'Position',[SV(1) Z(1) SV(2) Z(end)],'FaceColor',[1 1 0 0.4],'EdgeColor','k')
    plot(ax{11},v1_0,Zmdl,'sb','MarkerFaceColor','b');
    drawnow
end
% [HA0,r0,v1,ax] = CLASS_ClassicPointInv.Invert_v1(HA0,r0,v1_0,Zmdl,Zinv,Z,dZ,ao,ObsDta,options1,ax,[14]);
options1.StepTolerance = 1e-4;
NumLegendreCoeff = 30;
[HA0,r0,v1,ax] = CLASS_Legendre.Invert_Legendre_v1...
                              (NumLegendreCoeff,HA0,r0,v1_0,Zmdl,Zinv,Z,dZ,ao,ObsDta,options1,ax,[14]);
v2 = v1+90; v2(v2>180) = v2(v2>180)-180;
plot(ax{11},v2,Zmdl,'.r','LineWidth',2,'MarkerSize',14),drawnow
%% -------------------------------- INVERT HA
if FigVis == 1
    cla(ax{10}) 
    plot(ax{10},EV(:,5),EV(:,1),'-sk','MarkerFaceColor','w');
    plot(ax{10},HA0,Zmdl,'sb','MarkerFaceColor','b'); 
    drawnow
end
options1.StepTolerance = 1;
NumLegendreCoeff = 10;
[HA00,r0,v1,ax] = CLASS_Legendre.Invert_Legendre_HA...
                                   (NumLegendreCoeff,HA0,r0,v1,Zmdl,Zinv,Z,dZ,ao,ObsDta,options1,ax,[14]);
%% -------------------------------- INVERT r
if FigVis == 1
   cla(ax{12})
    plot(ax{12},r0,Zmdl,'sb','MarkerFaceColor','b'); 
    drawnow
end
options1.StepTolerance = 1;
NumLegendreCoeff = 10;
% [HA00,r00,v1,ax] = CLASS_ClassicPointInv.Invert_r(HA00,r0,v1,Zmdl,Zinv,Z,dZ,ao,ObsDta,options1,ax,[5]);
[HA00,r00,v1,ax] = CLASS_Legendre.Invert_Legendre_r...
                                 (NumLegendreCoeff,HA00,r0,v1,Zmdl,Zinv,Z,dZ,ao,ObsDta,options1,ax,[5]);
%% -------------------------------- GET EIGEN VALUES
if FigVis == 1
    cla(ax{9})
    plot(ax{9},EV(:,2),EV(:,1),'-sk','MarkerFaceColor','y');
    plot(ax{9},EV(:,3),EV(:,1),'-sk','MarkerFaceColor','c');
    plot(ax{9},EV(:,4),EV(:,1),'-sk','MarkerFaceColor','m');
    cla(ax{10})
    plot(ax{10},EV(:,5),EV(:,1),'-sk','MarkerFaceColor','w');
    plot(ax{10},HA0,Zmdl,'sb','MarkerFaceColor','b');
    plot(ax{10},HA00,Zmdl,'sg','MarkerFaceColor','g');
    cla(ax{12})
    plot(ax{12},r0,Zmdl,'sb','MarkerFaceColor','b'); 
    plot(ax{12},r00,Zmdl,'sg','MarkerFaceColor','g'); 
    drawnow
end
[EigVal,HA,r] = FUNC_GetEigenValues(ax,HA00,r00,Zmdl,ResZmdl,EV);
%% -------------------------------- OPTIMIZED MODEL
for i=5:8
    cla(ax{i});
end
tic
OptimizedModel = CLASS_FM.BeginForwardModel(Zmdl,HA,r,v1,dZ,0);
toc
EstPar = OptimizedModel.Dta;
estPA_HH = EstPar{5};
estPA_HV = EstPar{7};
estCP = EstPar{14};
estdPsi = EstPar{18};
ax = CLASS_InvPlot.UpdateInversionPlot(ax,estPA_HH,estPA_HV,estCP,estdPsi,Z,[],ao,[],[],FigVis);
fg.InvertHardcopy = 'off';
%% -------------------------------- Save parameters and plot
InversionOutput.MetaData= MetaData;
InversionOutput.OP_ObserevedData = RadarDta;
InversionOutput.OP_OptimizedModel = OptimizedModel;
InversionOutput.Bed = Bed;
InversionOutput.coreEV = EV;
InversionOutput.AxOut = AxOut;
InversionOutput.HA0 = HA0;
InversionOutput.r0 = r0;
InversionOutput.v1_0 = v1_0;
InversionOutput.v1 = v1;
InversionOutput.v2 = v2;
InversionOutput.r00 = r00;
InversionOutput.HA00 = HA00;
InversionOutput.EigVal = EigVal;
InversionOutput.HA = HA;
InversionOutput.r = r;
InversionOutput.Zmdl = Zmdl;
InversionOutput.Ztop = Ztop;
save(LabelName+"_InversionResults.mat",'InversionOutput');
print(fg,LabelName+"_InversionFigure.png",'-dpng','-r300');
end
%%
function [MHB] = funcHBDC(i,maxRange)
Ax_hb_All{1} = [166]; %1-W18 (not set)
Ax_hb_All{2} = [136]; %2-W12 (not set)
Ax_hb_All{3} = [173]; %3-W6 (set) (not set)
Ax_hb_All{4} = [185]; %4-W4.5 (not set)
Ax_hb_All{5} = []; %5-W3.5 (strange)
Ax_hb_All{6} = [162]; %6-W2.5 (not set) 
Ax_hb_All{7} = [176]; %7-W1.5 (not set)
Ax_hb_All{8} = [177]; %8-W1  (not set)
Ax_hb_All{9} = [124]; %9-W0.5 (not set)
Ax_hb_All{10} = [115]; %10-E0 (not set)
Ax_hb_All{11} = [198]; %11-EPICA (not set)
Ax_hb_All{12} = [183]; %12-E0.5 (not set)
Ax_hb_All{13} = [160]; %13-E1 (not set)
Ax_hb_All{14} = [174]; %14-E1.5 (not set)
Ax_hb_All{15} = [198]; %15-E2 (not set)
Ax_hb_All{16} = [195]; %16-E3 (not set)
Ax_hb_All{17} = [138]; %17-E4.5 (not set)
Ax_hb_All{18} = [135]; %18-E6 (not set)
Ax_hb_All{19} = [165]; %19-E9 (not set)
Ax_hb_All{20} = [194]; %20-E12 (not set)
Ax_hb_All{21} = [177]; %21-E18 (set)
MHB = Ax_hb_All{i};
MHB = [MHB' ; maxRange];
end
