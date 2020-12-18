clear
close all
clc
[DataPath,RadarData,InvPath,ps] = FUNC_ApRES_PathFix;
tm = string(datetime(now,'ConvertFrom','datenum'));
%% -------------------------------- Data Path
% Path of ApRES data
mainDataPath = '/home/reza/TIR_Brussels/';
% mainDataPath = '/mnt/esd01/docs/mershadi/WorkDirectory/ApRES/Data/RadarData_ApRES/TIR_Brussels/';
% mainDataPath = '\\134.2.5.43\esd01\docs\mershadi\WorkDirectory\ApRES\Data\RadarData_ApRES\TIR_Brussels\';
%% -------------------------------- Input Parameters
nm = 'p14';
DtaDir = strcat(mainDataPath,'2018',nm);
%% -------------------------------- Input Parameters
maxRange = 601;
BedRange = [500 600];
dA = 1;
ao = 0:dA:179; 
Rot_ao = 0;
C_DepthWin = maxRange * 0.01;
C_ConvWin = maxRange * 0.1;
DenoisingFlag.PA = ["1", "MovingAverage" , string(maxRange*0.05) ;
                                      "0", "Conv1D" , string(maxRange*0.1) ;
                                      "0", "Conv2D" , string(maxRange*0.05) ;
                                      "0", "DenoisePCA" , string(1)];
DenoisingFlag.PD = ["1", "MovingAverage" , string(maxRange*0.05) ;
                                      "0", "Conv1D" , string(maxRange*0.01) ;
                                      "0", "Conv2D" , string(maxRange*0.01) ;
                                      "0", "DenoisePCA" , string(1)];

[hh,vh,hv,vv,Z,dZ,maz,f,Bed,MetaData] = FUNC_ReadApRES(DtaDir,maxRange,BedRange);
%% Measurement correction
% temphh = hh;
% tempvh = vh;
% temphv= hv;
% tempvv = vv;
% hh = tempvv;
% vh = temphv;
% hv = tempvh;
% vv = temphh;
%%
GeoRef = 90; % Polarization plane correction
[~,izmx] = min(abs(Z-maxRange));
Z = Z(1:izmx,:); 
hh = hh(1:izmx,:); 
vh = hv(1:izmx,:); 
hv = hv(1:izmx,:); 
vv = vv(1:izmx,:);

[HH,VV,HV,VH] = CLASS_S2P.AzimuthSynthesizer(hh,vh,hv,vv,ao,GeoRef);
Dta = CLASS_S2P.Signal2Param(HH,VH,HV,VV,Z,ao,f,C_DepthWin,C_ConvWin,DenoisingFlag,"radar");

HHPA = Dta{5};
HVPA = Dta{7};
Phi = Dta{14};
Psi = Dta{18};

pltdim = [0.1,0.1,0.4,0.6];
w2plt = [5 14 7 18];
[fg,ax,cb] = CLASS_FixedPlot.AdvancePlot([],Dta,ao,Z,w2plt,Bed,pltdim);
for i = 2:length(ax)
yticklabels(ax{i},'')
ylabel(ax{i},[])
end

fg.InvertHardcopy = 'off';
print(fg,strcat(nm,'.png'),'-dpng','-r300')


% cb = colorbar(ax{1});
% cb.Label.String = '\deltaP_{HH} [-]';
% cb = colorbar(ax{2});
% cb.Label.String = '\phi_{HHVV} [rad]';
% cb = colorbar(ax{3});
% cb.Label.String = '\deltaP_{HV} [-]';
% cb = colorbar(ax{4});
% cb.Label.String = '\Psi [-]';
% 
% yticklabels(ax{2},'')
% ylabel(ax{2},[])
% yticklabels(ax{3},'')
% ylabel(ax{3},[])
% yticklabels(ax{4},'')
% ylabel(ax{4},[])
% 
% set(ax{1},'FontSize',20)
% set(ax{2},'FontSize',20)
% set(ax{3},'FontSize',20)
% set(ax{4},'FontSize',20)
% title(ax{1},['\fontsize{14}HH Power anomaly'])
% title(ax{2},['\fontsize{14}HHVV Coherence phase'])
% title(ax{3},['\fontsize{14}HV Power anomaly'])
% title(ax{4},['\fontsize{14}Scaled phase derivative'])
%%
% figure,plot(20*log10(abs(hh)),Z),set(gca,'YDIR','reverse')