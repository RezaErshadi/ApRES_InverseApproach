clear;
close all;
clc;
[DataPath,RadarData,sp] = FUNC_ApRES_PathFix;
% ---------------------------------------
Sensitive_param = "v1";
SL = 1;
% ---------------------------------------
dZ = 1;
hb = [500 1000]';
v1 = [0 0]';
HA = [0.2 0.2]';
r = [1 1]';
s2nr =0;
pltdim = [0.1,0.1,0.8,0.7];
w2plt = [5 7 14 18];
%%
fg = [];
i = 1;
while true
    switch Sensitive_param
        case "v1"
            P = 0:5:180;
            v1(SL) = P(i);
        case "HA"
            P = 0:0.05:0.5;
            HA(SL) = P(i);
        case "r"
            P = -30:5:30;
            r(SL) = P(i);
    end
    OP = CLASS_FM.BeginForwardModel(hb,HA,r,v1,dZ,s2nr);
    Dta = OP.Dta;
    ao = OP.ao;
    Z = OP.Z;
    [fg,ax,cb] = CLASS_FixedPlot.AdvancePlot(fg,Dta,ao,Z,w2plt,[],pltdim);
    drawnow
    fm(i) = getframe(fg);
    if i == length(P)
        break
    end
    i = i+1;
end
FUNC_SaveVideo(fm,2)