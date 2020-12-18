clear;
close all;
clc;
[dataPath] = FUNC_ApRES_StartLine; % fix the path for data and functions
sv = 0;
%% Linear transmitted signal
Ex0 = pi;
Ey0 = pi;
t = 0;
z=0:0.1:50; % depth : direction of propagation

[ExT,EyT] = FUNC_emw(0,Ex0,Ey0,t,z);
phi = -2*pi:pi/16:2*pi;
Ey0 = -2*pi:pi/16:2*pi;

CLASS_FixedPlot.SetFigureSize(0,0.1,0.5,0.8);
fg = gcf;
for i = 1:length(phi)
plot3(ExT,EyT,z,'b');
hold on
[ExR,EyR] = FUNC_emw(0,Ex0,Ey0(i),t,z);
plot3(ExR,EyR,z,'r');
hold off
set(gca,'ZDir','reverse');
grid on
xlabel('Amplitude')
ylabel('Amplitude')
zlabel('Depth')
ylim([-2*pi 2*pi])
xlim([-2*pi 2*pi])
zlim([0 max(z)])
% view(90,90)
% title("Phase shift: "+string(phi(i)/pi)+" \pi",'color','k')
set(gca,'FontSize',16,'Layer','top');
drawnow
pause(1)
F(i) = getframe(fg);
end
% FUNC_SaveVideo(F,1);


% FUNC_SetFigureSize(0,0.1,0.5,0.8);
% fg = gcf;
% for i = 1:length(phi)
% plot3(ExT,EyT,z,'b');
% hold on
% [ExR,EyR] = FUNC_emw(phi(i),Ex0,Ey0,t,z);
% plot3(ExR,EyR,z,'r');
% ExRT = ExR+ExT;
% EyRT = EyR+EyT;
% plot3(ExRT,EyRT,z,'k');
% hold off
% set(gca,'ZDir','reverse');
% grid on
% % xlabel('Amplitude')
% ylabel('Amplitude')
% zlabel('Depth')
% ylim([-2*pi 2*pi])
% % xlim([-pi pi])
% zlim([0 max(z)])
% view(90,0)
% title("Phase shift: "+string(phi(i)/pi)+" \pi",'color','k')
% set(gca,'FontSize',16,'Layer','top');
% drawnow
% pause(1)
% F(i) = getframe(fg);
% end
% FUNC_SaveVideo(F,1);






function [ExT,EyT] = FUNC_emw(phi,Ex0,Ey0,t,z)
    f = 300e6;
    w = 2 * pi * f; %(2*pi/T) angular frequency
    lambda = 5;
    k = 2*pi/lambda; % magnitude of the propagation vector
    wt = w*t;
    ExT=Ex0.*cos( (wt) - (k.*z) );
    EyT=Ey0.*cos( (wt) - (k.*z)+phi );
end