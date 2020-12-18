clear
close all
clc
[DataPath,RadarData,InvPath,ps]= FUNC_ApRES_PathFix;
%% Data path
mainDataPath = strcat(DataPath,ps,'RadarData_ApRES/EDML_AWI/PpRES_Cores/');
DtaDir = strcat(mainDataPath,'PpRES_LD01');
%%
maxRange = 2001;
dA = 22.5;
ao = 0:dA:179; 
C_DepthWin = maxRange * 0.01;
C_ConvWin = maxRange * 0.01;
DenoisingFlag.PA = ["1", "MovingAverage" , string(maxRange*0.05) ;
                                      "0", "Conv1D" , string(maxRange*0.01) ;
                                      "2", "Conv2D" , string(maxRange*0.05) ;
                                      "0", "DenoisePCA" , string(1)];
DenoisingFlag.PD = ["1", "MovingAverage" , string(maxRange*0.05) ;
                                      "0", "Conv1D" , string(maxRange*0.01) ;
                                      "0", "Conv2D" , string(maxRange*0.01) ;
                                      "0", "DenoisePCA" , string(1)];
%%
[hh,vh,hv,vv,Z,dZ,maz,f,Bed,MetaData] = FUNC_ReadApRES(DtaDir,maxRange,[]);
GeoRef = str2double(MetaData(MetaData(:,1) == "TR2TN",2));
GeoRef = GeoRef+ 90;
%% Measurement correction
% temphh = hh;
% tempvv = vv;
% hh = tempvv;
% vv = temphh;
 %% 
 % apply corrections for clockwise  (instead of counter clockwise)
 % use the whole HH to get VV

HHc = [hh(:,1) flip(hh(:,2:end),2)];
VHc = [hv(:,1) flip(hv(:,2:end),2)];
HVc = [hv(:,1) flip(hv(:,2:end),2)];
VVc = -[vv(:,1) flip(vv(:,2:end),2)];

aa = 1;
HH0 = HHc(:,aa); 
VH0 = VHc(:,aa); 
HV0 = HVc(:,aa); 
VV0 = VVc(:,aa); 

nn = 2;
HHc = [HHc(:,nn:end) HHc(:,1:nn-1)];
VHc = [VHc(:,nn:end) VHc(:,1:nn-1)];
HVc = [HVc(:,nn:end) HVc(:,1:nn-1)];
VVc = [VVc(:,nn:end) VVc(:,1:nn-1)];

pltdim = [0.1,0.1,0.4,0.6];
w2plt = [5 7 14];
%%
Dta = CLASS_S2P.Signal2Param(HHc,VHc,HVc,VVc,Z,maz,f,C_DepthWin,C_ConvWin,DenoisingFlag,"radar");
HHPAm = Dta{5};
HVPAm = Dta{7};
Phim = Dta{14};
Psim = Dta{18};

[fg,ax,cb] = CLASS_FixedPlot.AdvancePlot([],Dta,ao,Z,w2plt,Bed,pltdim);
for i = 2:length(ax)
yticklabels(ax{i},'')
ylabel(ax{i},[])
end

cb = colorbar(ax{1});
cb.Label.String = 'P_{HH} [-]';
cb = colorbar(ax{2});
cb.Label.String = 'P_{HV} [-]';
cb = colorbar(ax{3});
cb.Label.String = '\phi_{HHVV} [rad]';
% cb = colorbar(ax{4});
% cb.Label.String = '\Psi [-]';
yticklabels(ax{2},'')
ylabel(ax{2},[])
yticklabels(ax{3},'')
ylabel(ax{3},[])
% yticklabels(ax{4},'')
% ylabel(ax{4},[])
set(ax{1},'FontSize',14)
set(ax{2},'FontSize',14)
set(ax{3},'FontSize',14)
% set(ax{4},'FontSize',20)
title(ax{1},['\fontsize{14}(a1) HH Power anomaly'])
title(ax{2},['\fontsize{14}(a2) HV Power anomaly'])
title(ax{3},['\fontsize{14}(a3) HHVV Coherence phase'])
% title(ax{4},['\fontsize{14}Scaled phase derivative'])


% xticks(ax{1},[0 45 90 135])
% xticks(ax{2},[0 45 90 135])
% xticks(ax{3},[0 45 90 135])

%%
% ao = 0:1:179;
% [HH,VV,HV,VH] = CLASS_S2P.AzimuthSynthesizer(HH0,VH0,HV0,VV0,ao,GeoRef);
% Dta = CLASS_S2P.Signal2Param(HH,VH,HV,VV,Z,ao,f,C_DepthWin,C_ConvWin,DenoisingFlag,"radar");
% HHPAs = Dta{5};
% HVPAs = Dta{7};
% Phis = Dta{14};
% Psis = Dta{18};
% [fg,ax,cb] = CLASS_FixedPlot.AdvancePlot([],Dta,ao,Z,w2plt,Bed,pltdim);
% for i = 2:length(ax)
% yticklabels(ax{i},'')
% ylabel(ax{i},[])
% end
% 
% cb = colorbar(ax{1});
% cb.Label.String = 'P_{HH} [-]';
% cb = colorbar(ax{2});
% cb.Label.String = 'P_{HV} [-]';
% cb = colorbar(ax{3});
% cb.Label.String = '\phi_{HHVV} [rad]';
% % cb = colorbar(ax{4});
% % cb.Label.String = '\Psi [-]';
% yticklabels(ax{2},'')
% ylabel(ax{2},[])
% yticklabels(ax{3},'')
% ylabel(ax{3},[])
% % yticklabels(ax{4},'')
% % ylabel(ax{4},[])
% set(ax{1},'FontSize',14)
% set(ax{2},'FontSize',14)
% set(ax{3},'FontSize',14)
% % set(ax{4},'FontSize',20)
% title(ax{1},['\fontsize{14}(a1) HH Power anomaly'])
% title(ax{2},['\fontsize{14}(a2) HV Power anomaly'])
% title(ax{3},['\fontsize{14}(a3) HHVV Coherence phase'])
% % title(ax{4},['\fontsize{14}Scaled phase derivative'])
%%
fg.InvertHardcopy = 'off';

set(fg,'Units','inches');
screenposition = get(fg,'Position');
set(fg,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
print(fg,'LD01_measured225.pdf','-dpdf','-painters','-r300')
% print(fg,'LD01_synthesized225.pdf','-dpdf','-painters','-r300')
% print(fg,'LD01_synthesized1.pdf','-dpdf','-painters','-r300')
set(gcf,'Units','normalized');



% print(fg,"LD01_measured22.5.png",'-dpng','-r300')
% print(fg,"LD01_synthesized22.5.png",'-dpng','-r300')
% print(fg,"LD01_synthesized1.png",'-dpng','-r300')