clear
% close all
clc
[DataPath,RadarData,InvPath,ps] = FUNC_ApRES_PathFix;
tm = string(datetime(now,'ConvertFrom','datenum'));
%% -------------------------------- Data Path
% Path of ApRES data
mainDataPath = '/mnt/esd01/docs/mershadi/WorkDirectory/ApRES/Data/RadarData_ApRES/Concordia_BAS/';
% ff = ["W18","W12","W06","W4d5","W2d5","W1d5","W1d0","W0d5","E0",...
%         "EPICA",...
%         "E0d5","E1d0","E1d5","E02","E03","E4d5","E06","E09","E12","E18"];
% ss = ["1","2","3","4","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21"];
DtaDir = strcat(mainDataPath,'E03');
%% -------------------------------- Input Parameters
maxRange = 2001;
BedRange = [];
dA = 1;
ao = 0:dA:179; 
C_DepthWin = maxRange * 0.05;
C_ConvWin = maxRange * 0.05;
DenoisingFlag.PA = ["1", "MovingAverage" , string(maxRange*0.05) ;
                                      "0", "Conv1D" , string(maxRange*0.1) ;
                                      "2", "Conv2D" , string(maxRange*0.05) ;
                                      "0", "DenoisePCA" , string(1)];
DenoisingFlag.PD = ["1", "MovingAverage" , string(maxRange*0.05) ;
                                      "0", "Conv1D" , string(maxRange*0.01) ;
                                      "0", "Conv2D" , string(maxRange*0.01) ;
                                      "0", "DenoisePCA" , string(1)];

[hh,vh,hv,vv,Z,dZ,maz,f,Bed,MetaData] = FUNC_ReadApRES(DtaDir,maxRange,BedRange);
%% Measurement correction
% temphh = hh;
% tempvv = vv;
% hh = tempvv;
% vv = temphh;
%%
[~,izmx] = min(abs(Z-maxRange));
Z = Z(1:izmx,:); 

hh = hh(1:izmx,:); 
vh = hv(1:izmx,:); 
hv = hv(1:izmx,:); 
vv = vv(1:izmx,:);

% GeoRef = 90; % Polarization plane correction
GeoRef = str2double(MetaData(MetaData(:,1) == "TR2TN",2));
GeoRef = GeoRef +90;

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
% fg.InvertHardcopy = 'off';
% print(fg,"EPICA.png",'-dpng','-r300')


