clear all
close all
clc
%%
load('E1E2E3_Schmidt.mat')
% E2E1 = (0:0.001:0.5);
% E3E2 = (0:0.001:1);
% [X,Y] = meshgrid(E2E1,E3E2);
% for i = 1:length(E2E1)
%     tic
%     parfor j = 1:length(E3E2)
%         x = E2E1(i);
%         y = E3E2(j);
%         e3 = (1+(2*(y))+x)/3;
%         e2 = e3-y;
%         e1 = e2-x;
%         check1 = (e1 < 0 || e2 < 0 || e3 < 0);
%         check2 = (e1 > 1 || e2 > 1 || e3 > 1);
%         check3 = (e3<e1);
%         check4 = (e3<e2);
%         check5 = (e2<e1);
%         sumcheck = sum([check1,check2,check3,check4,check5]);
%         if sumcheck ~= 0
%             e1 = nan;
%             e2 = nan;
%             e3 = nan;
%         end
%         E1(j,i) = e1;
%         E2(j,i) = e2;
%         E3(j,i) = e3;
%     end
%     toc
% end
%%
cm = getPyPlot_cMap('seismic',21);
figure,
fg = gcf;
set(fg,'position',[1922         822        1918         694])
ax1 = subplot(1,12,1:3);
s = pcolor(X,Y,E1);
set(s, 'EdgeColor', 'none');
xlabel('E2-E1')
ylabel('E3-E2')
cb1 = colorbar;
cb1.Ticks = [linspace(min(E1(:)),max(E1(:)),3)];
cb1.TickLabels = {string(round(cb1.Ticks(1),2))+" (min)",string(round(cb1.Ticks(2),2)),string(round(cb1.Ticks(3),2))+" (max)"};
colormap(cm)
caxis([0 1])
title('E1')
set(ax1,'FontSize',16)
pbaspect([1 2 1])
hold on
plot(0,0,'ok','MarkerSize',20,'MarkerFaceColor',[0.5 0.5 0.5]);
plot(0,1,'ok','MarkerSize',20,'MarkerFaceColor',[0.5 0.5 0.5]);
plot(0.5,0,'ok','MarkerSize',20,'MarkerFaceColor',[0.5 0.5 0.5]);

ax2 = subplot(1,12,5:7);
s = pcolor(X,Y,E2);
set(s, 'EdgeColor', 'none');
xlabel('E2-E1')
% ylabel('E3-E2')
cb2 = colorbar;
cb2.Ticks = [linspace(min(E2(:)),max(E2(:)),3)];
cb2.TickLabels = {string(round(cb2.Ticks(1),2))+" (min)",string(round(cb2.Ticks(2),2)),string(round(cb2.Ticks(3),2))+" (max)"};
colormap(cm)
caxis([0 1])
title('E2')
set(ax2,'FontSize',16)
pbaspect([1 2 1])

ax3 = subplot(1,12,9:11);
s = pcolor(X,Y,E3);
set(s, 'EdgeColor', 'none');
xlabel('E2-E1')
% ylabel('E3-E2')
cb3 = colorbar;
cb3.Ticks = [linspace(min(E3(:)),max(E3(:)),3)];
cb3.TickLabels = {string(round(cb3.Ticks(1),2))+" (min)",string(round(cb3.Ticks(2),2)),string(round(cb3.Ticks(3),2))+" (max)"};
colormap(cm)
caxis([0 1])
title('E3')
set(ax3,'FontSize',16)
pbaspect([1 2 1])
%%
A(1,:) = X(1,:);
A(2,:) = Y(1,:);
A(3,:) = E1(1,:);
A(4,:) = E2(1,:);
C(1,:) = X(:,1)';
C(2,:) = flip(Y(:,1)');
C(3,:) = flip(E1(:,1)');
C(4,:) = flip(E2(:,1)');
for i = 1:size(Y,1)
    t = E1(i,:);
    t(isnan(t)) = [];
    j2 = length(t);
    B(1,i) = X(i,j2);
    B(2,i) = Y(i,1);
    B(3,i) = E1(i,j2);
    B(4,i) = E2(i,j2);
end
TT = [A,B,C];
TT(3:4,:) = rescale(TT(3:4,:),0.02,0.15);
ti = round(linspace(1,size(TT,2),500),0);
T(1,:) = TT(1,ti);
T(2,:) = TT(2,ti);
T(3,:) = TT(3,ti);
T(4,:) = TT(4,ti);
%%
cnt = [0.32,0.82];
circ = viscircles(ax1,[cnt(1) cnt(2)],0.15);
for i = 1:size(T,2)
    p = plot(ax1,T(1,i),T(2,i),'or','MarkerSize',10,'MarkerFaceColor','r');
    [xlp,ylp] =ellipse(T(3,i),T(4,i),0,cnt(1),cnt(2));
    [sx,isx] = sort(xlp);
    sy = ylp(isx);
    isy1 = find(sy>=cnt(2));
    sx1 = sx(isy1);
    Usy1 = sy(isy1);
    Dsy1 = Usy1-cnt(2);
    Lsy1 = cnt(2)-Dsy1;
    for j = 1:length(sx1)
        r(:,j) = Lsy1(j) + (Usy1(j)-Lsy1(j))*rand(3,1);
    end
    rp = plot(ax1,sx1,r,'.b');
    fab = plot(ax1,xlp,ylp,'-k','linewidth',2);
    F(i) = getframe(fg);
    drawnow
    delete(p)
    delete(fab)
    delete(rp)
end
% FUNC_SaveVideo(F,20);
%%
delete(circ)
delete(fab)
delete(rp)
delete(p)









