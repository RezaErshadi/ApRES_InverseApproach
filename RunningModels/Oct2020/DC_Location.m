clear;
close all;
clc;
[DataPath,RadarData,InvPath,ps] = FUNC_ApRES_PathFix;
%% ---------------------------------------
FilePath_DCLocation = strcat(DataPath,ps,'SiteInfo',ps,'DC_AntennaInfo.csv');
T = readtable(FilePath_DCLocation);
%%
TxLong = [table2array(T(:,5)) table2array(T(:,6)) table2array(T(:,7))];
TxLat = [table2array(T(:,9)) table2array(T(:,10)) table2array(T(:,11))];
RxLong = [table2array(T(:,13)) table2array(T(:,14)) table2array(T(:,15))];
RxLat = [table2array(T(:,17)) table2array(T(:,18)) table2array(T(:,19))];

TxLongDD = TxLong(:,1) + TxLong(:,2)./60 + TxLong(:,3)/3600;
TxLatDD = TxLat(:,1) + TxLat(:,2)./60 + TxLat(:,3)/3600;
RxLongDD = RxLong(:,1) + RxLong(:,2)./60 + RxLong(:,3)/3600;
RxLatDD = RxLat(:,1) + RxLat(:,2)./60 + RxLat(:,3)/3600;

[TxRxDistance,Az] = distance(TxLatDD,TxLongDD,RxLatDD,RxLongDD,wgs84Ellipsoid);
TxRxDistance = round(TxRxDistance,1);
Az = round(Az,1);
Az(Az>180) = Az(Az>180) -180;

T(:,20) = array2table(TxLongDD);
T(:,21) = array2table(TxLatDD);
T(:,22) = array2table(RxLongDD);
T(:,23) = array2table(RxLatDD);
T(:,24) = array2table(TxRxDistance);
T(:,25) = array2table(Az);
%%
FilePath_DCLocation2 = strcat(DataPath,ps,'SiteInfo',ps,'DC_AntennaInfo2.csv');
writetable(T,FilePath_DCLocation2);
