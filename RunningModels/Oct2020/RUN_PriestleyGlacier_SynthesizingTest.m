clear
close all
clc
[DataPath,RadarData,InvPath,ps] = FUNC_ApRES_PathFix; % Fix all the neccessary paths
tm = string(datetime(now,'ConvertFrom','datenum')); % Time of the code
%% Data path
loc = 7;
mainPath = '/mnt/esd01/docs/mershadi/WorkDirectory/ApRES/Data/RadarData_ApRES/PriestleyGlacier_Otago/';
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
        % possible combination:
        % hh -vv vh vh
        % -hh vv hv hv
    case 6
        DtaLabel = 'Site1_BAS_20min';
        % possible combination:
        % hh -vv vh vh
        % -hh vv hv hv
    case 7
        DtaLabel = 'Site2_BAS_7min'; 
        % possible combination:
        % hh vv hv hv
        % hh vv -vh -vh
    case 8
        DtaLabel = 'Site2_BAS_20min';
        % possible combination:
        % hh vv hv hv
        % hh vv -vh -vh
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
nm = ["hh" "vv" "hv" "vh"];
pltdim = [0.1,0.1,0.5,0.7];
CLASS_FixedPlot.SetFigureSize(pltdim(1),pltdim(2),pltdim(3),pltdim(4));    fg = gcf;   fg.Color = "white";
for i = 1:4
    eval('SignalPower = 20*log10(abs('+nm(i)+'));')
    eval('SignalPhase = angle('+nm(i)+');')
    phsuw = rescale(unwrap(SignalPhase),-1,1);
    subplot(2,4,i)
    plot(SignalPower,Z)
    hold on
    plot([min(SignalPower) max(SignalPower)],[Bed Bed],'--r','linewidth',2)
    set(gca,'Ydir','reverse')
    ylim([Z(1) Z(end)])
    xlim([min(SignalPower) max(SignalPower)])
    xlabel('Power [dB]')
    ylabel('Depth [m]')
    grid on
    title("Power "+nm(i))
    subplot(2,4,i+4)
    plot(SignalPhase,Z)
    hold on
    plot(phsuw,Z,'-g')
    plot([min(SignalPhase) max(SignalPhase)],[Bed Bed],'--r','linewidth',2)
    set(gca,'Ydir','reverse')
    ylim([Z(1) Z(end)])
    xlim([min(SignalPhase) max(SignalPhase)])
    grid on
    title("Phase "+nm(i))
    xlabel('Phase [rad]')
    ylabel('Depth [m]')
end
print(fg,strcat('Site2_BAS_7min_Power&Phase'),'-dpng','-r300');
%%
% trshld = 250;
% [hh,axhh] = func_Signal2Noise(hh,Z,maxRange,Bed,trshld,1);
% [vh,axvh] = func_Signal2Noise(vh,Z,maxRange,Bed,trshld,1);
% [hv,axhv] = func_Signal2Noise(hv,Z,maxRange,Bed,trshld,1);
% [vv,axvv] = func_Signal2Noise(vv,Z,maxRange,Bed,trshld,1);
%%
if contains(DtaDir,'BAS') % Synthesizing
    pltdim = [0.1,0.1,0.7,0.4];
    w2plt = [5 7 14 18];
% replace
% load('SIte1_Comb(Az-BAS).mat');
% hh = vv_az90;
% vv = vv_az0;
%----
    HH0 = -hh(:,m_az==0); 
    VH0 = vh(:,m_az==0);
    HV0 = vh(:,m_az==0); 
    VV0 = vv(:,m_az==0); 
    [HH,VV,HV,VH] = CLASS_S2P.AzimuthSynthesizer(HH0,VH0,HV0,VV0,ao,0);
    Dta = CLASS_S2P.Signal2Param(HH,VH,HV,VV,Z,ao,f,C_DepthWin,C_ConvWin,DenoisingFlag,"radar");
    [fg,ax,cb] = CLASS_FixedPlot.AdvancePlot([],Dta,ao,Z,w2plt,Bed,pltdim);
%     fctr = combvec([1 -1],[1 -1],[1 -1],[1 -1]);
%     for i = 1:size(fctr,2)
%         hh0 = fctr(1,i)*hh(:,m_az==0); 
%         vh0 = fctr(2,i)*vh(:,m_az==0);
%         hv0 = fctr(3,i)*hv(:,m_az==0); 
%         vv0 = fctr(4,i)*vv(:,m_az==0); 
%         [HH,VV,HV,VH] = CLASS_S2P.AzimuthSynthesizer(hh0,vh0,hv0,vv0,ao,0);
%         Dta = CLASS_S2P.Signal2Param(HH,VH,HV,VV,Z,ao,f,C_DepthWin,C_ConvWin,DenoisingFlag,"radar");
%         [fg3,ax3,cb3] = CLASS_FixedPlot.AdvancePlot([],Dta,ao,Z,w2plt,[],pltdim);
%         print(fg3,strcat('Site2_BAS_20min_Synthesized_',string(i),'(hh,vh,hv,vv)'),'-dpng','-r300');
%     end
%     for i = 1:size(fctr,2)
%         hh0 = fctr(1,i)*hh(:,m_az==0); 
%         vh0 = fctr(2,i)*hv(:,m_az==0);
%         hv0 = fctr(3,i)*hv(:,m_az==0); 
%         vv0 = fctr(4,i)*vv(:,m_az==0); 
%         [HH,VV,HV,VH] = CLASS_S2P.AzimuthSynthesizer(hh0,vh0,hv0,vv0,ao,0);
%         Dta = CLASS_S2P.Signal2Param(HH,VH,HV,VV,Z,ao,f,C_DepthWin,C_ConvWin,DenoisingFlag,"radar");
%         [fg3,ax3,cb3] = CLASS_FixedPlot.AdvancePlot([],Dta,ao,Z,w2plt,[],pltdim);
%         print(fg3,strcat('Site2_BAS_20min_Synthesized_',string(i),'(hh,hv,hv,vv)'),'-dpng','-r300');
%     end
%     for i = 1:size(fctr,2)
%         hh0 = fctr(1,i)*hh(:,m_az==0); 
%         vh0 = fctr(2,i)*vh(:,m_az==0);
%         hv0 = fctr(3,i)*vh(:,m_az==0); 
%         vv0 = fctr(4,i)*vv(:,m_az==0); 
%         [HH,VV,HV,VH] = CLASS_S2P.AzimuthSynthesizer(hh0,vh0,hv0,vv0,ao,0);
%         Dta = CLASS_S2P.Signal2Param(HH,VH,HV,VV,Z,ao,f,C_DepthWin,C_ConvWin,DenoisingFlag,"radar");
%         [fg3,ax3,cb3] = CLASS_FixedPlot.AdvancePlot([],Dta,ao,Z,w2plt,[],pltdim);
%         print(fg3,strcat('Site2_BAS_20min_Synthesized_',string(i),'(hh,vh,vh,vv)'),'-dpng','-r300');
%     end
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
print(fg,FigName,'-dpng','-r300');

%%
function [specCor,ax] = func_Signal2Noise(specCor,range,maxRange,BedDepth,trshld,plt)
    SignaldB = abs(specCor);
    NoiseFloor=max(SignaldB(range>maxRange-trshld));
    NoiseFloordb=20*log10(NoiseFloor);
    ax = [];
    if plt == 1
        figure; hold all;
        ax = gca;
        plot(20*log10(abs(specCor)),range,'.r');
        set(ax,'Ydir','reverse')
        set(ax,'Color',[0.6 0.6 0.6]);
        if ~isempty(BedDepth)
            plot([NoiseFloordb,NoiseFloordb],[0,BedDepth],'-b','linewidth',2);
            plot([NoiseFloordb,NoiseFloordb],[BedDepth maxRange],'--b','linewidth',2);
            plot([min(20*log10(SignaldB)),max(20*log10(SignaldB))],[BedDepth,BedDepth],'--k','linewidth',2);
        else
            plot([NoiseFloordb,NoiseFloordb],[0,maxRange],'--k','linewidth',2);
        end
        plot(20*log10(abs(specCor(abs(specCor)>NoiseFloor))),range(abs(specCor)>NoiseFloor),'.g');
        title('Noise | Signal')
        xlabel('Power [dB]')
        ylabel('Depth [m]')
        grid on
        xlim([min(20*log10(SignaldB)),max(20*log10(SignaldB))])
        ylim([range(1) range(end)])
    end
    specCor(SignaldB<NoiseFloor)=nan;
end