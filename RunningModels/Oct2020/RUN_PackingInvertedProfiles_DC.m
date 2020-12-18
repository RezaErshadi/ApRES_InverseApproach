clear
close all
clc
[DataPath,RadarData,InvPath,ps] = FUNC_ApRES_PathFix;
tm = string(datetime(now,'ConvertFrom','datenum'));
%%
DataFolder = strcat('Corrected_HHVV',ps,'DC_FinalInversion',ps,'Data');
fp = strcat(InvPath,ps,DataFolder);
DataList = dir(string(fp)+ps+"*.mat");
for i = 1:length(DataList)
    cn(i,1) = string(DataList(i).name);
end
idel_EPICA = find(contains(cn,'EPICA'));
idel_MAT = find(contains(cn,'INV.mat'));
idel = [idel_EPICA ; idel_MAT];
clear('cn');
DataList(idel,:) = [];
for i = 1:length(DataList)
    cn = DataList(i).name;
    dsh = strfind(cn,'_');
    sID(1,i) = str2double(cn(1:dsh(1)-1));
end
[sID,isID] = sort(sID,'descend');
DataList = DataList(isID);
%%
tic
for i = 1:length(DataList)
    tic
    fprintf(string(i)+'\n')
dtapath = strcat(DataList(i).folder,ps,DataList(i).name);
load(dtapath);
%%
%MethaData
MD(:,i) = InversionOutput.MetaData(:,2);
SiteID(:,i) = MD(1:2,i);
SiteElevation(1,i) = str2double(MD(11,i));
SiteBed(:,i) = str2double(MD(12:13,i));
% observed coherence phase
if i == 1
    Z = InversionOutput.OP_ObserevedData.Z;
    dZ = InversionOutput.OP_ObserevedData.dZ;
    ao = InversionOutput.OP_ObserevedData.ao;
    Zmdl = InversionOutput.Zmdl;
end
CP(:,:,i) = InversionOutput.OP_ObserevedData.Dta{14};
% Inverted parameters
v1(:,i) = InversionOutput.v1;
v2(:,i) = InversionOutput.v2;
HA(:,i) = InversionOutput.HA;
r(:,i) = InversionOutput.r;
lambda1(:,i) = InversionOutput.EigVal(:,1);
lambda2(:,i) = InversionOutput.EigVal(:,2);
lambda3(:,i) = InversionOutput.EigVal(:,3);
    toc
end
toc
save(strcat(fp,ps,'DC_INV.mat'),...
    'MD','SiteID','SiteElevation','SiteBed','ao','Z','dZ','CP','Zmdl','v1','v2','HA','r','lambda1','lambda2','lambda3')







