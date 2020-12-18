clear
close all
clc
[DataPath,RadarData,ps] = FUNC_ApRES_PathFix;
%%
ForwardModelFlag = 0;
%% Data path
loc = 3;
switch loc
    case 1
%         DtaDir = ...
%         '/mnt/esd01/docs/mershadi/WorkDirectory/ApRES/Data/RadarData_ApRES/Concordia_BAS/E0';  %10 E0 DomeC
        DtaDir = ...
        '/mnt/esd01/docs/mershadi/WorkDirectory/ApRES/Data/RadarData_ApRES/Concordia_BAS/EPICA';  %11 EPICA DomeC
    case 2
        DtaDir = ...
        '/mnt/esd01/docs/mershadi/WorkDirectory/ApRES/Data/RadarData_ApRES/EDML_AWI/PpRES_Cores/PpRES_LD01';  %14 LD01 EDML
    case 3
        DtaDir = ...
        '/mnt/esd01/docs/mershadi/WorkDirectory/ApRES/Data/RadarData_ApRES/PriestleyGlacier_Otago/Site_1_20mins'; %Priestley
end
%%
maxRange = 1001;
dA = 1;
ao = 0:dA:179; 
Rot_ao = 0;
C_DepthWin = maxRange * 0.025;
C_ConvWin = maxRange * 0.025;
DenoisingFlag.PA = ["1", "MovingAverage" , string(maxRange*0.025) ;
                                      "0", "Conv1D" , string(maxRange*0.01) ;
                                      "1", "Conv2D" , string(maxRange*0.025) ;
                                      "0", "DenoisePCA" , string(1)];
DenoisingFlag.PD = ["1", "MovingAverage" , string(maxRange*0.01) ;
                                      "0", "Conv1D" , string(maxRange*0.01) ;
                                      "0", "Conv2D" , string(maxRange*0.01) ;
                                      "0", "DenoisePCA" , string(1)];
pltdim = [0.1,0.1,0.9,0.4];
w2plt = [5 8 7 6 14];
%%
[hh,vv,hv,vh,Z,dZ,maz,f,Bed] = FUNC_ReadApRES(DtaDir,maxRange,[ ]);
%%
HH0 = hh(:,maz==0); 
VH0 = hv(:,maz==0);
HV0 = -hv(:,maz==0); 
VV0 = -vv(:,maz==0); 
[HH,VV,HV,VH] = CLASS_S2P.AzimuthSynthesizer(HH0,VH0,HV0,VV0,ao,0);
Dta = CLASS_S2P.Signal2Param(HH,VH,HV,VV,Z,ao,f,C_DepthWin,C_ConvWin,DenoisingFlag,"radar");
[fg1,ax1,cb1] = CLASS_FixedPlot.AdvancePlot([],Dta,ao,Z,w2plt,Bed,pltdim);
print(fg1,"Priestley_Site1_(HH,-VV,HV,-HV,StrongSymmetry).png",'-dpng','-r300');


% 
% [HH,VV,HV,VH] = CLASS_S2P.SignalDeveloper(HH0,VV0,HV0,HV0,ao+Rot_ao);
% dtaParams = CLASS_S2P.Signal2Param(HH,VV,HV,VH,Z,ao,f,NoiseWin,CohWin,ConvWin,dns,"r");
% [fg2,ax2,cb2] = CLASS_FixedPlot.AdvancePlot([],dtaParams,ao,Z,w2plt,[],pltdim);
% 
% [HH,VV,HV,VH] = CLASS_S2P.SignalDeveloper(HH0,VV0,VH0,VH0,ao+Rot_ao);
% dtaParams = CLASS_S2P.Signal2Param(HH,VV,HV,VH,Z,ao,f,NoiseWin,CohWin,ConvWin,dns,"r");
% [fg3,ax3,cb3] = CLASS_FixedPlot.AdvancePlot([],dtaParams,ao,Z,w2plt,[],pltdim);
