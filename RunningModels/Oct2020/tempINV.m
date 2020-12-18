clear
close all
clc
%%
load('DC_Flowdirection.mat')

X = DC(:,2);
Y = DC(:,1);
Az = deg2rad(DC(:,4));
Mag = DC(:,3);
[U,V] = pol2cart(Az,Mag);
figure,
compass(U,V)

load('EDML_Flowdirection.mat')
X = EDML(:,2);
Y = EDML(:,1);
Az = deg2rad(EDML(:,4));
Mag = EDML(:,3);
[U,V] = pol2cart(Az,Mag);
figure,
compass(U,V)