clear
close all
clc
[DataPath,RadarData,InvPath,ps] = FUNC_ApRES_PathFix;
tm = string(datetime(now,'ConvertFrom','datenum'));
%% Read surface and profile data (REMA & BEDMACHIB)
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
DataFolder = strcat('Corrected_HHVV',ps,'DC_FinalInversion',ps,'Data');
fp = strcat(InvPath,ps,DataFolder);
load(string(fp)+ps+"DC_INV.mat")
nData = length(SiteID);
%%
for i = 1:nData
temp1= char(SiteID(1,i));
SiteID(3,i) = temp1(2:end);
end
siteDist = str2double(SiteID(3,:));
cnt = find(siteDist == 0);
siteDist(1:cnt) = -siteDist(1:cnt);
tcklbl =  ["-18 (EAST)" "-12" "-9" "-6" "-4.5" "-3" "-2" "" "-1" "" ...
                      "0" "" "1" "" "2.5" "4.5" "6" "12" "18 (WEST)"];
%%
xx = min(siteDist):0.5:max(siteDist);
XX = [];
for i = 1:length(xx)
    XX = [XX repmat(xx(i),1,length(ao))]; 
end
aa = 1:length(XX);
Zobs = (3300:-dZ:(3300-2201))';
Zinv = (3300:-mean(diff(Zmdl)):(3300-2201))';
matCP = nan(length(Zobs),length(XX));
matHA = nan(length(Zinv),length(XX));
matV2 = nan(length(Zinv),length(XX));
%%
for i = 1:nData
    a = find(XX==siteDist(i));
    xcnt(i) = a(length(a)/2);
    tempCP = CP(:,:,i);
    tempHA = repmat(lambda2(:,i)-lambda1(:,i),1,size(tempCP,2));
    tempV2 = repmat(v2(:,i),1,size(tempCP,2));
    [~,j11] = min(abs(SiteElevation(i)-Zobs)); % top obs
    [~,j22] = min(abs(SiteElevation(i)-Zinv)); % top inv
    for j = 1:length(a)
        matCP(j11:j11+size(tempCP,1)-1,a(j)) = tempCP(:,j);
        matHA(j22:j22+size(tempHA,1)-1,a(j)) = tempHA(:,j); 
        matV2(j22:j22+size(tempV2,1)-1,a(j)) = tempV2(:,j);
    end
    tempCP = matCP(:,a(1):a(end));
    [~,m1] = min(abs(pi-tempCP(:)));
    Elvmat = repmat(Zobs,1,length(a));
    NodeElv(1,i) = Elvmat(m1);
    disp(i)
end
%%
Bprof = [BP(:,1) nan(length(BP),1) BP(:,4)];
for i = 1:length(SiteBed)
    [~,t1(1,i)] = min(abs(SiteBed(1,i)-Bprof(:,3)));
end
t1(2,:) = xcnt;
for i = 1:length(t1)
    Bprof(t1(1,i),2) = t1(2,i);
end
for i = 1:length(t1)-1
Bprof(t1(1,i):t1(1,i+1),2) = rescale(Bprof(t1(1,i):t1(1,i+1),1),Bprof(t1(1,i),2) ,Bprof(t1(1,i+1),2));
end
Bprof(isnan(Bprof(:,2)),:) = [];
%%
bkgrndCM = [0.65 0.65 0.65];
caxV2 = [0 180];
caxHA = [0 0.5];

cm1 = getPyPlot_cMap('seismic',100);
ll = -0.01; % girdle strength lower limit
cmlength = length(-0.5:abs(ll):0.5);
cm21 = getPyPlot_cMap('cubehelix',50);
cm22 = repmat([1 1 1],200,1);
cm2 = [cm21 ; cm22];
cm3 = flip(getPyPlot_cMap('cubehelix',90));

pltdim = [0.5000    0.0256    0.5003    0.9213];
CLASS_FixedPlot.SetFigureSize(pltdim(1),pltdim(2),pltdim(3),pltdim(4));    
fg1 = gcf;   fg1.Color = "white";
ax1 = subplot(10,1,1);
ax2 = subplot(10,1,[2 4]);
ax3 = subplot(10,1,[5 7]);
ax4 = subplot(10,1,[8 10]);
%  --------------------------------
xSprof = linspace(aa(90),aa(13050),length(SP));
yyaxis(ax1,'right')
plot(ax1,Bprof(:,2),Bprof(:,3),'-','MarkerSize',10,'linewidth',2)
hold(ax1,'on')
plot(ax1,xcnt,SiteBed(2,:),'xr','MarkerSize',10,'linewidth',3)
ylim(ax1,[-200 200])
ylabel(ax1,'Bed [m asl]')
yyaxis(ax1,'left')
plot(ax1,xSprof',SP(:,4),'.','MarkerSize',10,'linewidth',3)
% hold(ax1,'on')
% plot(ax1,xcnt,SiteElevation(1,:),'x','MarkerSize',10,'linewidth',3)
ylim(ax1,[3220 3240])
ylabel(ax1,'Surface [m asl]')
xticks(ax1,xcnt)
xticklabels(ax1,[])
ax1.XGrid = 'on';
xlim(ax1,[aa(1) aa(end)])
% set(ax1,'FontSize',12)
%  --------------------------------
% yyaxis(ax1,'right')
% plot(ax1,xcnt,SiteBed(2,:),'x-','MarkerSize',10,'linewidth',3)
% ylim(ax1,[-200 200])
% ylabel(ax1,'Bed')
% yyaxis(ax1,'left')
% plot(ax1,xcnt,SiteElevation(1,:),'x-','MarkerSize',10,'linewidth',3)
% ylim(ax1,[3220 3240])
% ylabel(ax1,'Surface')
% xticks(ax1,xcnt)
% xticklabels(ax1,[])
% ax1.XGrid = 'on';
% xlim(ax1,[aa(1) aa(end)])
% set(ax1,'FontSize',12)
% --------------------------------
imagesc(ax2,aa,-Zobs,matCP)
hold(ax2,'on')
plot(ax2,xcnt,-NodeElv,':g','linewidth',3)
set(ax2,'YDIR','reverse')
colormap(ax2,[bkgrndCM; cm1])
cb2 = colorbar(ax2,'Location','eastoutside');
cb2.Label.String = "\phi_{HHVV} [rad]";
cb2.Ticks = [-3.14 -2 -1 0 1 2 3.14];
cb2.TickLabels = {'-\pi','-2','-1','0','1','2','\pi'};
caxis(ax2,[-pi pi])
xticks(ax2,xcnt)
xticklabels(ax2,[])
grid(ax2,'on')
ylabel(ax2,'Elevation [m asl]')
yticklabels(ax2,-double(string(cell2mat(get(ax2,'YTickLabel')))))
% --------------------------------
s = pcolor(ax3,aa,Zinv,matHA);
shading(ax3,'flat');
s.EdgeColor = 'none';
ax3.Color = bkgrndCM;
colormap(ax3,cm2)
cb3 = colorbar(ax3,'Location','eastoutside');
cb3.Label.String = "\Delta\lambda [-]";
caxis(ax3,caxHA)
xticks(ax3,xcnt)
xticklabels(ax3,[])
grid(ax3,'on')
ylabel(ax3,'Elevation [m asl]')
% --------------------------------
s = pcolor(ax4,aa,Zinv,matV2);
shading(ax4,'flat');
s.EdgeColor = 'none';
ax4.Color = bkgrndCM;
colormap(ax4,cm3)
cb4 = colorbar(ax4,'Location','eastoutside');
cb4.Label.String = "v2 [deg]";
cb4.Ticks = 0:45:180;
caxis(ax4,caxV2)
xticks(ax4,xcnt)
xticklabels(ax4,tcklbl)
grid(ax4,'on')
xlabel(ax4,"Distance to the Dome [km]");
ylabel(ax4,'Elevation [m asl]')
% --------------------------------
cbp2 = [0.910590277777776 0.61287988422576 0.0126091734628955 0.225781235600796];
cbp3 = get(cb3,'Position');
cbp4 = get(cb4,'Position');
set(cb2,'Position',[cbp2(1) cbp2(2) cbp3(3) cbp2(4)]);
set(cb3,'Position',[cbp2(1) cbp3(2) cbp3(3) cbp3(4)]);
set(cb4,'Position',[cbp2(1) cbp4(2) cbp4(3) cbp4(4)]);
set(ax1,'FontSize',14)
set(ax2,'FontSize',14)
set(ax3,'FontSize',14)
set(ax4,'FontSize',14)
fg1.InvertHardcopy = 'off';
%%
% set(fg1,'Units','inches');
% screenposition = get(fg1,'Position');
% set(fg1,...
%         'PaperPosition',[0 0 screenposition(3:4)],...
%         'PaperSize',[screenposition(3:4)]);
% print(fg1,"FinalResult",'-dpdf','-painters','-r300')
% set(fg1,'Units','normalized');
% --------------------------------
% yax3 = double(string(cell2mat(get(ax3,'YTickLabel'))));
% yax4 = double(string(cell2mat(get(ax4,'YTickLabel'))));
% yticklabels(ax3,-yax3);
% yticklabels(ax4,-yax4);
% --------------------------------
print(fg1,"FinalResult.png",'-dpng','-r300')






