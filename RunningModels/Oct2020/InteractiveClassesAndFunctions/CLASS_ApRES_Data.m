classdef CLASS_ApRES_Data
    % All the related functions to load ApRES data
    methods(Static)
function Concordia = Read_Concordia(dataPath,Si,Zmx,Rot_ao,dA,NoiseWin,CohWin,ConvWin)
    fid = fopen('ApRES_All_Info.csv','r');
    c = textscan(fid,'%f %s %f %f %f %f %f %f %f','Delimiter',',','HeaderLines', 1);
    fclose(fid);
    ApRES_name = c{2};
    ApRES = [c{1} c{3} c{4} c{5} c{6} c{7} c{8} c{9}];
    nData = length(ApRES);
    %%
    dns = 1; % denoising
    ao = 0:dA:179;      aom=(ao(1:end-1)+ao(2:end))/2; % antenna orientation
    %% Reading the data
    DtaDir = strcat(dataPath,['/ApRES/Concordia']);
    [HH0,VV0,HV0,Z,dZ,f,DtaName] = FUNC_ReadRadarDataConcordia(Si,Zmx,DtaDir);
    D2 = -1.*(Z-ApRES(Si,6)); % set the ApRES top at the surface elevation
    %% Developing the signal
    if DtaName{1} == 11
        Rot_ao = 90;
    end
    [HH,VV,HV,VH] = CLASS_S2P.SignalDeveloper(HH0,VV0,HV0,HV0,ao+Rot_ao);
    dtaParams = CLASS_S2P.Signal2Param(HH,VV,HV,VH,Z,ao,f,NoiseWin,CohWin,ConvWin,dns,"r");
    %% outpur
    Concordia.f = f;
    Concordia.ao = ao;
    Concordia.Z = Z;
    Concordia.dZ = dZ;
    Concordia.Zmx = Zmx;
    Concordia.OrgHH = [HH];
    Concordia.OrgVV = [VV];
    Concordia.OrgHV = [HV];
    Concordia.OrgVH = [VH];
    Concordia.dtaParams = dtaParams;
end
%% -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
function EDML = Read_EDML(dataPath,Si,Zmx,Rot_ao,NoiseWin,CohWin,ConvWin)
    dns = 1; % denoising
    ao = 0:1:179;      aom=(ao(1:end-1)+ao(2:end))/2; % antenna orientation
    [HH0,VV0,HV0,Z,dZ,f] = FUNC_ReadRadarData(dataPath,3,Si,"y",Zmx);
    [HH,VV,HV,VH] = CLASS_S2P.SignalDeveloper(HH0,VV0,HV0,HV0,ao+Rot_ao);
    dtaParams = CLASS_S2P.Signal2Param(HH,VV,HV,VH,Z,ao,f,NoiseWin,CohWin,ConvWin,dns,"r");
    %% outpur
    EDML.f = f;
    EDML.ao = ao;
    EDML.Z = Z;
    EDML.dZ = dZ;
    EDML.Zmx = Zmx;
    EDML.OrgHH = [HH];
    EDML.OrgVV = [VV];
    EDML.OrgHV = [HV];
    EDML.OrgVH = [VH]; 
    EDML.dtaParams = dtaParams;     
end
%% -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
function [bd] = ApRES_Bed(HH0,HV0,VV0,Z,BedRange)
    bnHH = fmcw_findbed(Z(:,1),abs(HH0(:,1)),BedRange,'xcor');
    bnVV = fmcw_findbed(Z(:,1),abs(VV0(:,1)),BedRange,'xcor');
    bnHV = fmcw_findbed(Z(:,1),abs(HV0(:,1)),BedRange,'xcor');
    bd = mean([Z(bnHH) Z(bnVV) Z(bnHV)]);
end
%% -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
function EV = LoadEigenValues(FP)
    fid = fopen(FP,'r');
    Dta = textscan(fid, '%f %f %f %f %s', 'headerlines', 1, 'Delimiter',',');
    fclose(fid);
    ZEV = Dta{:,1}; %depth
    E1 = Dta{:,2}; %E1
    E2 = Dta{:,3}; %E2
    E3 = Dta{:,4}; %E2
    GirdleStrength = E2-E1;
    PoleStrength = E3-E2;
    EV = [ZEV,E1,E2,E3,GirdleStrength,PoleStrength];
end
%% -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*




%% -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*


%% -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    end
end