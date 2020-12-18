clear all
close all
clc
%%
load('E1E2E3_Schmidt.mat')
%%
cm = getPyPlot_cMap('seismic',21);
figure,
fg = gcf;
set(fg,'position',[1922         822        1918         694])
ax1 = subplot(1,3,1);
s = pcolor(ax1,X,Y,E1);
set(s, 'EdgeColor', 'none');
xlabel('\lambda2-\lambda1')
ylabel('\lambda3-\lambda2')
cb1 = colorbar;
cb1.Ticks = [min(E1(:)),max(E1(:))];
% cb1.TickLabels = {string(round(cb1.Ticks(1),2))+" (min E1)",string(round(cb1.Ticks(2),2)),string(round(cb1.Ticks(3),2))+" (max E1)"};
cb1.TickLabels = {string(round(cb1.Ticks(1),2))+" (min \lambda1)",string(round(cb1.Ticks(2),2))+" (max \lambda1)"};
colormap(cm)
caxis([0 1])
title('\lambda1')
set(ax1,'FontSize',16)
pbaspect([1 2 1])
hold on
% plot(0,0,'ok','MarkerSize',20,'MarkerFaceColor',[0.5 0.5 0.5]);
% plot(0,1,'ok','MarkerSize',20,'MarkerFaceColor',[0.5 0.5 0.5]);
% plot(0.5,0,'ok','MarkerSize',20,'MarkerFaceColor',[0.5 0.5 0.5]);
xticks(ax1,[0:0.1:0.5]);
yticks(ax1,[0:0.1:1]);
set(ax1,'color',[0.5 0.5 0.5])
grid(ax1,'on')
set(ax1,'layer','top')

ax2 = subplot(1,3,2);
s = pcolor(ax2,X,Y,E2);
hold on
set(s, 'EdgeColor', 'none');
xlabel('\lambda2-\lambda1')
% ylabel('E3-E2')
cb2 = colorbar;
cb2.Ticks = [min(E2(:)),max(E2(:))];
cb2.TickLabels = {string(round(cb2.Ticks(1),2))+" (min \lambda2)",string(round(cb2.Ticks(2),2))+" (max \lambda2)"};
colormap(cm)
caxis([0 1])
title('\lambda2')
set(ax2,'FontSize',16)
xticks(ax2,[0:0.1:0.5]);
yticks(ax2,[0:0.2:1]);
yticklabels(ax2,[]);
set(ax2,'color',[0.5 0.5 0.5])
grid(ax2,'on')
pbaspect([1 2 1])
set(ax2,'layer','top')

ax3 = subplot(1,3,3);
s = pcolor(ax3,X,Y,E3);
hold on
set(s, 'EdgeColor', 'none');
xlabel('\lambda2-\lambda1')
% ylabel('E3-E2')
cb3 = colorbar;
cb3.Ticks = [0 0.33 0.5 1];
% cb3.Label.String = "Eigenvalue";
cb3.Ticks = [min(E3(:)),max(E3(:))];
cb3.TickLabels = {string(round(cb3.Ticks(1),2))+" (min \lambda3)",string(round(cb3.Ticks(2),2))+" (max \lambda3)"};
colormap(cm)
caxis([0 1])
title('\lambda3')
set(ax3,'FontSize',16)
xticks(ax3,[0:0.1:0.5]);
yticks(ax3,[0:0.2:1]);
yticklabels(ax3,[]);
set(ax3,'color',[0.5 0.5 0.5])
grid(ax3,'on')
pbaspect([1 2 1])
set(ax3,'layer','top')
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
    p1 = plot(ax1,T(1,i),T(2,i),'or','MarkerSize',15,'MarkerFaceColor','g');
    p2 = plot(ax2,T(1,i),T(2,i),'or','MarkerSize',15,'MarkerFaceColor','g');
    p3 = plot(ax3,T(1,i),T(2,i),'or','MarkerSize',15,'MarkerFaceColor','g');
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
    fg.InvertHardcopy = 'off';
    F(i) = getframe(fg);
    drawnow
    delete(p1)
    delete(p2)
    delete(p3)
    delete(fab)
    delete(rp)
end
% FUNC_SaveVideo(F,20);
%%
delete(circ)
delete(fab)
delete(rp)
delete(p1)
delete(p2)
delete(p3)
%%
% fg.InvertHardcopy = 'off';
% print(fg,"FabricType.png",'-dpng','-r300')
FUNC_SaveVideo(F,20);






