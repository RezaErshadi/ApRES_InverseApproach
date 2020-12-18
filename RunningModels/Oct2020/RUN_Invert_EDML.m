clear
close all
clc
[DataPath,RadarData,InvPath,ps] = FUNC_ApRES_PathFix;
tm = string(datetime(now,'ConvertFrom','datenum'));
%% -------------------------------- Data Path
% Path of ApRES data
DtaDir = '/home/reza/OneDrive(bwedu)/WorkDirectory/ApRES/Data/RadarData_ApRES/EDML_AWI/PpRES_Cores/PpRES_LD01';  %14 LD01 EDML
% Path of EigenValues data
EvDir = 'Data/EigenValues/EDML_EigenValues_Weikusat.txt';
Invert = 1;
%% -------------------------------- Input Parameters
maxRange = 2001;
BedRange = [];
ResZmdl = 50; % Model Depth Resolution
dA = 1;
ao = 0:dA:179; 
Rot_ao = 0;
C_DepthWin = maxRange * 0.01;
C_ConvWin = maxRange * 0.01;
DenoisingFlag.PA = ["1", "MovingAverage" , string(maxRange*0.05) ;
                                      "0", "Conv1D" , string(maxRange*0.01) ;
                                      "2", "Conv2D" , string(maxRange*0.01) ;
                                      "0", "DenoisePCA" , string(1)];
DenoisingFlag.PD = ["1", "MovingAverage" , string(maxRange*0.05) ;
                                      "0", "Conv1D" , string(maxRange*0.01) ;
                                      "0", "Conv2D" , string(maxRange*0.01) ;
                                      "0", "DenoisePCA" , string(1)];
pltdim = [0.1,0.1,0.4,0.6];
w2plt = [5 7 9 14 17 18];
FigVis = 1;
%% -------------------------------- Read the Data
[hh,vh,hv,vv,Z,dZ,maz,f,Bed,MetaData] = FUNC_ReadApRES(DtaDir,maxRange,BedRange);
%% Apply orientation correction
hh = [hh(:,1) flip(hh(:,2:end),2)];
vh = [hv(:,1) flip(hv(:,2:end),2)];
hv = [hv(:,1) flip(hv(:,2:end),2)];
vv = -[vv(:,1) flip(vv(:,2:end),2)];
[~,izmx] = min(abs(Z-maxRange));
Z = Z(1:izmx,:); hh = hh(1:izmx,:); vh = vh(1:izmx,:); hv = hv(1:izmx,:); vv = vv(1:izmx,:);
SiteNumber = str2double(MetaData(MetaData(:,1) == "SiteNumber",2));
GeoRef = str2double(MetaData(MetaData(:,1) == "TR2TN",2));
GeoRef = GeoRef+ 90;
%%
HH0 = hh(:,maz==0); 
VH0 = hv(:,maz==0);
HV0 = hv(:,maz==0); 
VV0 = vv(:,maz==0); 
[HH,VV,HV,VH] = CLASS_S2P.AzimuthSynthesizer(HH0,VH0,HV0,VV0,ao,GeoRef);
Dta = CLASS_S2P.Signal2Param(HH,VH,HV,VV,Z,ao,f,C_DepthWin,C_ConvWin,DenoisingFlag,"radar");
RadarDta.f = f;
RadarDta.ao = ao;
RadarDta.Z = Z;
RadarDta.dZ = dZ;
RadarDta.Zmx = max(Z);
RadarDta.Dta = Dta;
FoldData = split(DtaDir,ps);
LabelName = strcat(InvPath,ps,FoldData{end-1},'_',FoldData{end},'_',tm);
%% --------------------------------
% function [InversionOutput,fg,ax] = InvertApRES(EvDir,SiteIndex,RadarDta,ResZmdl,LabelName,FigVis)
ao = RadarDta.ao;
Z = RadarDta.Z;
dZ = RadarDta.dZ;
Zmx = RadarDta.Zmx;
ObsDta = RadarDta.Dta;
PAHH = ObsDta{5};
PAHV = ObsDta{7};
CP = ObsDta{14};
Psi = ObsDta{18};
EV = CLASS_ApRES_Data.LoadEigenValues(EvDir);
MHB = funcHBDC(SiteNumber,Zmx);
Zmdl = (ResZmdl:ResZmdl:Zmx)'; Zmdl(end) = Zmx;
%%
cp1 = CP;
aomat = repmat(ao,length(Z),1);
Zmat = repmat(Z,1,length(ao));
cp2 = find(abs(pi-cp1)<1e-3);
NodeInfo(1,:) = aomat(cp2);
NodeInfo(2,:) = Zmat(cp2);
cp3 = cp1(1:end-1,:).*cp1(2:end,:);
cp3 = [cp1(1,:) ; cp3];
cp4 = find(cp3<0);
BorderInfo(1,:) = aomat(cp4);
BorderInfo(2,:) = Zmat(cp4);
cp5 = zeros(size(CP));
cp5(cp4) = 1;
cp5 = [Z sum(cp5,2)];
cp5(cp5(:,1)<150 | cp5(:,1)>300,:) = [];
[~,imhb] = max(cp5(:,2));
Ztop = cp5(imhb,1);
%% --------------------------------
[fg,ax] = CLASS_InvPlot.InversionPlot(PAHH,PAHV,CP,Psi,Z,ao,[],FigVis);
% plot ice core Eigenvalues
plot(ax{9},EV(:,2),EV(:,1),'-sk','MarkerFaceColor','y');
plot(ax{9},EV(:,3),EV(:,1),'-sk','MarkerFaceColor','c');
plot(ax{9},EV(:,4),EV(:,1),'-sk','MarkerFaceColor','m');
% plot ice core lambda2-lambda1
plot(ax{10},EV(:,5),EV(:,1),'-sk','MarkerFaceColor','w');
AxOut = CLASS_Initials.PossiblePrincipalAxis(PAHV,ao);
% plot min of HV power anomaly on: dPHV and Psi
% plot(ax{2},BorderInfo(1,:),BorderInfo(2,:),'.k','linewidth',2)
plot(ax{3},AxOut.AxMin(:,1),Z,'.c','MarkerFaceColor','c')
plot(ax{3},AxOut.AxMin(:,2),Z,'.c','MarkerFaceColor','c')
plot(ax{4},AxOut.AxMin(:,1),Z,'.c','MarkerFaceColor','c')
plot(ax{4},AxOut.AxMin(:,2),Z,'.c','MarkerFaceColor','c')
%--------------------------------- Add North and surface velocity direction
% SV = [20 20];
% rectangle(ax{11},'Position',[SV(1) Z(1) SV(2) Z(end)],'FaceColor','y','EdgeColor','k')
% --------------------------------
j1 = 1;
for i = 1:size(MHB,1)
    [~,j2] = min(abs(MHB(i,1)-Zmdl));
    [~,ii] = min(abs(MHB(i,2)-AxOut.AxMean));
    v1_0(j1:j2,1) =  AxOut.AxMean(ii);
    j1 = j2+1;
end
% Plot initial v1
plot(ax{11},v1_0,Zmdl,'sb','MarkerFaceColor','b'); drawnow
% --------------------------------
HA0 = CLASS_Initials.PossibledE(Z,Zmdl,AxOut,Psi); % Initial Horizontal Anisotropy
% Plot initial HA
plot(ax{10},HA0,Zmdl,'sb','MarkerFaceColor','b'); drawnow 
% --------------------------------
r0 = zeros(size(Zmdl)); % Initial Reflection Ratio
% Plot initial r
plot(ax{12},r0,Zmdl,'sb','MarkerFaceColor','b'); drawnow 
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
%% -------------------------------- INVERT r
if FigVis == 1
   cla(ax{12})
    plot(ax{12},r0,Zmdl,'sb','MarkerFaceColor','b'); 
end
% options1.StepTolerance = 1e-4;
[HA0,r00,v1_0,ax] = CLASS_ClassicPointInv.Invert_r(HA0,r0,v1_0,Zmdl,Zinv,Z,dZ,ao,ObsDta,options1,ax,[5,7]);
% NumLegendreCoeff = 30;
% [HA0,r00,v1_0,ax] = CLASS_Legendre.Invert_Legendre_r(NumLegendreCoeff,HA0,r0,v1_0,Zmdl,Zinv,Z,dZ,ao,ObsDta,options1,ax,[5,7]);
%% -------------------------------- INVERT Ax
if FigVis == 1
    cla(ax{11})
%     rectangle(ax{11},'Position',[SV(1) Z(1) SV(2) Z(end)],'FaceColor',[1 1 0 0.4],'EdgeColor','k')
    plot(ax{11},v1_0,Zmdl,'sb','MarkerFaceColor','b');
    drawnow
end
options1.StepTolerance = 1e-4;
[HA0,r00,v1,ax] = CLASS_ClassicPointInv.Invert_v1(HA0,r00,v1_0,Zmdl,Zinv,Z,dZ,ao,ObsDta,options1,ax,[5,7]);
% NumLegendreCoeff = 30;
% [HA0,r00,v1,ax] = CLASS_Legendre.Invert_Legendre_v1(NumLegendreCoeff,HA0,r00,v1_0,Zmdl,Zinv,Z,dZ,ao,ObsDta,options1,ax,[5,7]);
v2 = v1+90; v2(v2>180) = v2(v2>180)-180;
plot(ax{11},v2,Zmdl,'.r','LineWidth',2,'MarkerSize',14); drawnow
%% -------------------------------- INVERT HA
if FigVis == 1
    cla(ax{10}) 
    plot(ax{10},EV(:,5),EV(:,1),'-sk','MarkerFaceColor','w');
    plot(ax{10},HA0,Zmdl,'sb','MarkerFaceColor','b'); 
    drawnow
end
HA00 = HA0;
% NumLegendreCoeff = 10;
% [HA00,r00,v1,ax] = CLASS_Legendre.Invert_Legendre_HA...
%                                     (NumLegendreCoeff,HA0,r00,v1,Zmdl,Zinv,Z,dZ,ao,ObsDta,options1,ax,[5,14]);
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
% SaveData = input("Should I save the results?\ny\nn\n---> ","s");
% if SaveData == "y"
    InversionOutput.MetaData= MetaData;
    InversionOutput.OP_ObserevedData = RadarDta;
    InversionOutput.OP_OptimizedModel = OptimizedModel;
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
    save(LabelName+"_InversionResults.mat",'InversionOutput');
    print(fg,LabelName+"_InversionFigure.png",'-dpng','-r300');
% end
%%
function [MHB] = funcHBDC(i,maxRange)
Ax_hb_All{14} = [[200;800;maxRange]  , [170;90;90]]; %1-W18 (not set)
MHB = Ax_hb_All{i};
end

