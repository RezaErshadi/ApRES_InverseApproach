clear
close all
clc
[DataPath,RadarData,InvPath,ps] = FUNC_ApRES_PathFix;
tm = string(datetime(now,'ConvertFrom','datenum'));
%%
load('E1E2E3_Schmidt.mat')
%%
cnt = [0 0];
distX = 0.2;
distY = 0.2;




xl = [0.22 0.28];
yl = [0.22 0.28];
xr = xl(1) + (xl(2)-xl(1)) .* rand(500,2);
yr = yl(1) + (yl(2)-yl(1)) .* rand(500,2);
plot(xr,yr,'.')
xlim([0 0.5])
ylim([0 0.5])





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
figure,
fg =gcf;
ax1 = gca;
cnt = [0,0];
% circ = viscircles(ax1,[cnt(1) cnt(2)],0.15);
for i = 1:size(T,2)
%     p = plot(ax1,T(1,i),T(2,i),'or','MarkerSize',10,'MarkerFaceColor','r');
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
%     F(i) = getframe(fg);
    xlim(ax1,[-0.5 0.5])
    ylim(ax1,[-0.5 0.5])
    drawnow
%     delete(p)
%     delete(fab)
%     delete(rp)
end