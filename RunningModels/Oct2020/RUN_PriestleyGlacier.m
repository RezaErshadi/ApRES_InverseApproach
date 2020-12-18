clear
close all
clc
[DataPath,RadarData,InvPath,ps] = FUNC_ApRES_PathFix; % Fix all the neccessary paths
tm = string(datetime(now,'ConvertFrom','datenum')); % Time of the code
%% Data path
loc = 1;
mainPath = '/mnt/esd01/docs/mershadi/WorkDirectory/ApRES/Data/RadarData_ApRES/PriestleyGlacier_Otago/';
% mainPath = '\\134.2.5.43\esd01\docs\mershadi\WorkDirectory\ApRES\Data\RadarData_ApRES\PriestleyGlacier_Otago\';
switch loc
    case 1
        DtaLabel = 'Site1_Azimuthal_0-170';
    case 2
        DtaLabel = 'Site1_Azimuthal_180-350';
    case 3
        DtaLabel = 'Site2_Azimuthal_0-170';
    case 4
        DtaLabel = 'Site2_Azimuthal_180-350'; 
    case 5
        DtaLabel = 'Site1_BAS_7min';
        fctr = [1 -1 -1 -1];
        % possible combination:
        % hh and vv opposite
        % hv and vh opposite
    case 6
        DtaLabel = 'Site1_BAS_20min';
        fctr = [1 -1 -1 -1];
    case 7
        DtaLabel = 'Site2_BAS_7min'; 
        fctr = [1 1 1 1];
        % possible combination:
        % hh and vv same
        % hv and vh opposite
    case 8
        DtaLabel = 'Site2_BAS_20min';
        fctr = [1 1 1 1];
    case 9
        DtaLabel = 'SiteAA_BAS_7min';
        fctr = [1 1 1 -1];
    case 10
        DtaLabel = 'SiteCC_BAS_7min';
        fctr = [1 1 1 -1];
end
DtaDir = strcat(mainPath,DtaLabel);
FigName = strcat(DtaLabel,'.png');
%%
maxRange = 1000;
BedRange = [700 maxRange];
dA = 1;
ao = 0:dA:179; 
Rot_ao = 0;
C_DepthWin = maxRange * 0.05;
C_ConvWin = maxRange * 0.05;
DenoisingFlag.PA = ["1", "MovingAverage" , string(maxRange*0.05) ;
                                      "0", "Conv1D" , string(maxRange*0.05) ;
                                      "0", "Conv2D" , string(maxRange*0.05) ;
                                      "0", "DenoisePCA" , string(1)];
DenoisingFlag.PD = ["0", "MovingAverage" , string(maxRange*0.01) ;
                                      "0", "Conv1D" , string(maxRange*0.01) ;
                                      "0", "Conv2D" , string(maxRange*0.01) ;
                                      "0", "DenoisePCA" , string(1)];
%%
[hh,vh,hv,vv,Z,dZ,m_az,f,Bed,MetaData] = FUNC_ReadApRES(DtaDir,maxRange,BedRange);
%%
% trshld = 250;
% [hh,axhh] = func_Signal2Noise(hh,Z,maxRange,Bed,trshld,1);
% [vh,axvh] = func_Signal2Noise(vh,Z,maxRange,Bed,trshld,1);
% [hv,axhv] = func_Signal2Noise(hv,Z,maxRange,Bed,trshld,1);
% [vv,axvv] = func_Signal2Noise(vv,Z,maxRange,Bed,trshld,1);
%%
if contains(DtaDir,'BAS') % Synthesizing
    pltdim = [0.1,0.1,0.7,0.4];
    w2plt = [5 14 7 18];
    HH0 = fctr(1)*hh(:,m_az==0); 
    VH0 = fctr(2)*hv(:,m_az==0);
    HV0 = fctr(3)*hv(:,m_az==0); 
    VV0 = fctr(4)*vv(:,m_az==0); 
    [HH,VV,HV,VH] = CLASS_S2P.AzimuthSynthesizer(HH0,VH0,HV0,VV0,ao,0);
    Dta = CLASS_S2P.Signal2Param(HH,VH,HV,VV,Z,ao,f,C_DepthWin,C_ConvWin,DenoisingFlag,"radar");
    [fg,ax,cb] = CLASS_FixedPlot.AdvancePlot([],Dta,ao,Z,w2plt,Bed,pltdim);
    caxis(ax{4},[0 1])
    cm = colormap(ax{4});
    cmlength = size(cm,1);
    cm2 = getPyPlot_cMap('cividis',cmlength);
    cm3 = getPyPlot_cMap('cool',cmlength);
    cm4 = [cm2 ; cm3];
    colormap(ax{4},cm4);
PAHV = Dta{7};
AxOut = CLASS_Initials.PossiblePrincipalAxis(PAHV,ao);
hold(ax{3},'on')
hold(ax{4},'on')
plot(ax{3},AxOut.AxMin(:,1),Z,'.g','MarkerFaceColor','g')
plot(ax{3},AxOut.AxMin(:,2),Z,'.g','MarkerFaceColor','g')
plot(ax{4},AxOut.AxMin(:,1),Z,'.g','MarkerFaceColor','g')
plot(ax{4},AxOut.AxMin(:,2),Z,'.g','MarkerFaceColor','g')
% print(fg,FigName,'-dpng','-r300');
elseif contains(DtaDir,'Azimuthal') % Azimuthal Method
    [PA_VV,PD_VV] = CLASS_S2P.AmpPhs(vv,m_az);
    PA.HH = nan;
    PA.VH = nan;
    PA.HV = nan;
    PD.HH = nan;
    PD.VH = nan;
    PD.HV = nan;
    PA.VV = PA_VV;
    PD.VV = PD_VV;
    [PA,PD] = CLASS_Denoising.RunDenoising(DenoisingFlag,PA,PD,dZ);
    PA_VV = PA.VV;
    PD_VV = PD.VV;
    nn = nan(size(PA_VV));
    for i = 1:18
        Dta{i} = nn;
    end
    Dta{8} = PA_VV;
    Dta{12} = PD_VV;
    pltdim = [0.1,0.1,0.15,0.5];
    w2plt = [8];
    [fg,ax,cb] = CLASS_FixedPlot.AdvancePlot([],Dta,m_az,Z,w2plt,[],pltdim);
end