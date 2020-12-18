clear;
close all;
clc;
[DataPath,RadarData,sp] = FUNC_ApRES_PathFix;
% use forward model
dZ = 1;
hb = [500 1000]';
Ax = [40 160 10]';
dE = [0.3 0.1 0.1]';
r = [3 -5 12]';
s2nr =0;
OP = CLASS_FM.BeginForwardModel(hb,dE,r,Ax,dZ,s2nr); 
%
f = OP.f;
ao = OP.ao;
Z = OP.Z;
dZ = OP.dZ;
Zmx = OP.Zmx;
Dta = OP.Dta;
if s2nr == 1
    DenoisingFlag=[]; C_DepthWin=dZ; C_ConvWin=dZ; 
else
    DenoisingFlag=[]; C_DepthWin=50; C_ConvWin=50; 
end
%
pltdim = [0.1,0.1,0.8,0.7];
w2plt = [5 7 14];
% [fg1,ax1,cb1] = CLASS_FixedPlot.AdvancePlot([],Dta,ao,Z,w2plt,[],pltdim);
%%
hh = Dta{1};
vh = Dta{2};
hv = Dta{3};
vv = Dta{4};
m_az = 0:10:170;
for i = 1:length(m_az)
    mHH(:,i) = hh(:,ao==m_az(i));
    mVH(:,i) = vh(:,ao==m_az(i));
    mHV(:,i) = hv(:,ao==m_az(i));
    mVV(:,i) = vv(:,ao==m_az(i));
end
m_Dta = CLASS_S2P.Signal2Param(mHH,mVH,mHV,mVV,Z,m_az,f,C_DepthWin,C_ConvWin,DenoisingFlag,"model");
if s2nr == 1
    m_Dta{5} = CLASS_Denoising.MovingAverage(m_Dta{5},50,dZ);
    m_Dta{7} = CLASS_Denoising.MovingAverage(m_Dta{7},50,dZ);
end
[fg1,ax1,cb1] = CLASS_FixedPlot.AdvancePlot([],m_Dta,m_az,Z,w2plt,[],pltdim);
print(fg1,'10deg_measured','-dpng','-r300');
%%
hh0 = mHH(:,m_az==0);
vh0 = mVH(:,m_az==0);
hv0 = mHV(:,m_az==0);
vv0 = mVV(:,m_az==0);
[HH,VV,HV,VH] = CLASS_S2P.AzimuthSynthesizer(hh0,vh0,hv0,vv0,ao,0);
Dta = CLASS_S2P.Signal2Param(HH,VH,HV,VV,Z,ao,f,C_DepthWin,C_ConvWin,DenoisingFlag,"model");
if s2nr == 1
    Dta{5} = CLASS_Denoising.MovingAverage(Dta{5},50,dZ);
    Dta{7} = CLASS_Denoising.MovingAverage(Dta{7},50,dZ);
end
[fg2,ax2,cb2] = CLASS_FixedPlot.AdvancePlot([],Dta,ao,Z,w2plt,[],pltdim);
print(fg2,'1deg_SynthesizedCorrect','-dpng','-r300');
%%
fctr = combvec([1 -1],[1 -1],[1 -1],[1 -1]);
for i = 1:size(fctr,2)
    hh0 = fctr(1,i)*mHH(:,m_az==0);
    vh0 = fctr(2,i)*mVH(:,m_az==0);
    hv0 = fctr(3,i)*mHV(:,m_az==0);
    vv0 = fctr(4,i)*mVV(:,m_az==0);
    [HH,VV,HV,VH] = CLASS_S2P.AzimuthSynthesizer(hh0,vh0,hv0,vv0,ao,0);
    Dta = CLASS_S2P.Signal2Param(HH,VH,HV,VV,Z,ao,f,C_DepthWin,C_ConvWin,DenoisingFlag,"model");
    [fg3,ax3,cb3] = CLASS_FixedPlot.AdvancePlot([],Dta,ao,Z,w2plt,[],pltdim);
    print(fg3,strcat('1deg_Synthesized_',string(i)),'-dpng','-r300');
end
%%
% mHH = dtaParams{1};
% mHV = dtaParams{2};
% mVV = dtaParams{3};
% mVH = dtaParams{4};
% %% sens test
% hh = dtaParams{1}(:,1);
% hv = dtaParams{2}(:,1);
% vv = dtaParams{3}(:,1);
% vh = dtaParams{4}(:,1);
% dns = 0; NoiseWin=dZ; CohWin=dZ; ConvWin=dZ; 
% fctr = combvec([1 -1],[1 -1]);
% for i = 1:size(fctr,2)
% HH0 = hh; 
% VV0 = vv; 
% HV0 = fctr(1,i)*hv; 
% VH0 = fctr(2,i)*vh; 
% 
% [HH,VV,HV,VH] = CLASS_S2P.SignalDeveloper(HH0,VV0,HV0,VH0,ao);
% dtaParams2 = CLASS_S2P.Signal2Param(HH,VV,HV,VH,Z,ao,f,NoiseWin,CohWin,ConvWin,dns,"s");
% [fg1,ax1,cb1] = CLASS_FixedPlot.AdvancePlot([],dtaParams2,ao,Z,w2plt,[],pltdim);
% drawnow
% end
% %%
% hh = dtaParams{1};
% hv = dtaParams{2};
% vv = dtaParams{3};
% vh = dtaParams{4};
% HH0 = hh(:,ao==0);
% VV0 = vv(:,ao==90);
% HV0 = hv(:,ao==0);
% VH0 = vh(:,ao==90);
% dns = 0; NoiseWin=dZ; CohWin=dZ; ConvWin=dZ; 
% [HH,VV,HV,VH] = CLASS_S2P.SignalDeveloper(HH0,VV0,HV0,VH0,ao);
% dtaParams2 = CLASS_S2P.Signal2Param(HH,VV,HV,VH,Z,ao,f,NoiseWin,CohWin,ConvWin,dns,"s");
% [fg1,ax1,cb1] = CLASS_FixedPlot.AdvancePlot([],dtaParams2,ao,Z,w2plt,[],pltdim);








