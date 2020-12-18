clear
close all
clc
[DataPath,RadarData,InvPath,ps] = FUNC_ApRES_PathFix;
tm = string(datetime(now,'ConvertFrom','datenum'));
%%
SurfaceFile  = strcat(DataPath,ps,'SiteInfo',ps,'DC_Profile_SurfaceElevation(REMA).tsv');
BedFile  = strcat(DataPath,ps,'SiteInfo',ps,'DC_Profile_BedElevation(BedMachin).tsv');
SurfaceData = tdfread(SurfaceFile);
SP(:,1) = SurfaceData.x00x2E0;
SP(:,2) = SurfaceData.x13416610x2E6011302285;
SP(:,3) = SurfaceData.x0x2D9017850x2E1423307826;
SP(:,4) = SurfaceData.x32240x2E860595703125;
BedData = tdfread(BedFile);
BP(:,1) = BedData.x00x2E0;
BP(:,2) = BedData.x13416610x2E6011302285;
BP(:,3) = BedData.x0x2D9017850x2E1423307826;
BP(:,4)= BedData.x160x2E996337890625;
%%