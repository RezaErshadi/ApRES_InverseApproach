clear;
close all;
clc;
[dataPath] = FUNC_ApRES_StartLine; % fix the path for data and functions
sv = 0;
%% Linear transmitted signal
phi = 0; % relative phase differencebetween the two components
Ex0 = pi;
Ey0 = pi;
f = 300e6;
w = 2 * pi * f; %(2*pi/T) angular frequency
t = 0:0.1:1;
lambda = 5;
k = 2*pi/lambda; % magnitude of the propagation vector
z=0:0.1:25; % depth : direction of propagation
wt = pi;

ExT=Ex0.*cos( (wt) - (k.*z) );
EyT=Ey0.*cos( (wt) - (k.*z)+phi );

%% Received signal
phi = -2*pi:pi/8:2*pi;
Ex0 = pi;
Ey0 = pi;
CLASS_FixedPlot.SetFigureSize(0,0.1,0.7,0.7);
fg = gcf;
fg.Color = 'black';
if min(Ey0) == max(Ey0)
    xylm = [-pi pi];
else
    xylm = [-pi pi];
end
zlm = [0 max(z)];
while true
    for n = 1:length(phi)
        ExR=Ex0.*cos( (wt) - (k.*z) );
        EyR=Ey0.*cos( (wt) - (k.*z)+phi(n) );
        %%
        subplot(2,2,1)
        plot3(ExT,EyT,z);
        set(gca,'ZDir','reverse');
        hold on
        plot3(ExR,EyR,z);
        set(gca,'ZDir','reverse');
        sbplt(xylm,zlm,135,0)
        hold off
        title("phase difference: "+string(phi(n)/pi)+" \pi",'color','y')

        subplot(2,2,2)
        plot3(ExT,EyT,z);
        set(gca,'ZDir','reverse');
        hold on
        plot3(ExR,EyR,z);
        set(gca,'ZDir','reverse');
        sbplt(xylm,zlm,0,90)
        hold off

        subplot(2,2,3)
        plot3(ExT,EyT,z);
        set(gca,'ZDir','reverse');
        hold on
        plot3(ExR,EyR,z);
        set(gca,'ZDir','reverse');
        sbplt(xylm,zlm,90,0)
        hold off
        
        subplot(2,2,4)
        plot3(ExT,EyT,z);
        set(gca,'ZDir','reverse');
        hold on
        plot3(ExR,EyR,z);
        set(gca,'ZDir','reverse');
        sbplt(xylm,zlm,30,30)
        hold off

        drawnow
        pause(0.5)
        
        F(n) = getframe(fg);
    end
    break
end
if sv == 1
    FUNC_SaveVideo(F,1);
end






%%
function []=sbplt(xylm,zlm,az,el)
    set(gca,'YDir','reverse');
    axis square;
    grid on;
    xlim(xylm)
    ylim(xylm)
    zlim(zlm)
    xlabel("X",'color','y')
    ylabel("Y",'color','y')
    zlabel("Depth",'color','y')
    set(gca,'XColor','y');
    set(gca,'YColor','c');
    set(gca,'ZColor','c');
    ax = gca;
    ax.GridColor = 'r';
    view(az, el);
end