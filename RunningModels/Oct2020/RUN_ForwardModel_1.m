clear;
close all;
clc;
[DataPath,RadarData,InvPath,ps] = FUNC_ApRES_PathFix;
tm = string(datetime(now,'ConvertFrom','datenum'));
% -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
figvis = 1;
dZ = 1;
s2nr = 0;
% [hb] = sort(CLASS_FM.RandBetween2Values(300,3000,[5,1]));
% [hb] = [500 1000];
% [dE] = [0.1 0.3]';
% [Ax] = [20 100]';
% [r] = [3 3 3 3 3]';
hb = (20:20:2000);
dE = nan(length(hb),1);
r = nan(length(hb),1);
Ax = nan(length(hb),1);

LC = CLASS_Legendre.FUNC_legendrefit(   [(0.33:-0.01:0.15)   (0.15:-0.01:0.12)]',5, 'inv');
E1 = CLASS_Legendre.FUNC_LegendreCoeff2Data([hb]',LC);

LC = CLASS_Legendre.FUNC_legendrefit(   [0.33:-0.01:0.13]',5, 'inv');
E2 = CLASS_Legendre.FUNC_LegendreCoeff2Data([hb]',LC);

% LC = CLASS_Legendre.FUNC_legendrefit(   [(0.1:0.01:0.2) (0.2:-0.01:0.1)]',5, 'inv');
% dE = CLASS_Legendre.FUNC_LegendreCoeff2Data([hb]',LC);

dE = E2-E1;
E3 = 1-E2-E1;

for i = 1:length(hb)-1
    E1T = E1(i);
    E2T = E2(i);
    E1B = E1(i+1);
    E2B = E2(i+1);
    r(i,1) = ((E2T-E2B).^2) ./  ((E1T-E1B).^2);
    r(i,1) = 20.*log10(r(i,1));
end
r(end) = r(end-1);

% LC = CLASS_Legendre.FUNC_legendrefit([(10:2:25) (25:-2:15)]',5, 'inv');
% Ax = CLASS_Legendre.FUNC_LegendreCoeff2Data([hb]',LC);
% Ax = repmat(25,length(hb),1);
[~,iax] = min(abs(hb-500)); 
Ax(1:iax) = 10;
Ax(iax+1:end) = 20;

% LC = CLASS_Legendre.FUNC_legendrefit(   [-5:1:15]',5, 'inv');
% r = CLASS_Legendre.FUNC_LegendreCoeff2Data([hb]',LC);

% r = [repmat(3,50,1) ; repmat(3,51,1)];
% [dE] = CLASS_FM.RandBetween2Values(0,0.5,[length(hb),1]);
% [Ax] = CLASS_FM.RandBetween2Values(0,179,[length(hb),1]);
% [r] = CLASS_FM.RandBetween2Values(-20,20,[length(hb),1]);
%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
OP_ObserevedData = CLASS_FM.BeginForwardModel(hb,dE,r,Ax,dZ,s2nr);  
ao = OP_ObserevedData.ao;
Z = OP_ObserevedData.Z;
dZ = OP_ObserevedData.dZ;
ObsDta = OP_ObserevedData.dtaParams;
PAHH = ObsDta{5};
PAHV = ObsDta{6};
PCP = ObsDta{16};
dEhor = ObsDta{19};
%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
% w2plt = [5,6,16,19];
% [fg,ax,cb] = CLASS_FixedPlot.AdvancePlot([],dtaParams,ao,Z,w2plt,[]); 
[fg,ax] = CLASS_InvPlot.InversionPlot(PAHH,PAHV,PCP,dEhor,Z,ao,[],figvis);

EV = [hb' E1 E2 E3 dE];
plot(ax{9},EV(:,2),EV(:,1),'sy','MarkerFaceColor','y');
plot(ax{9},EV(:,3),EV(:,1),'+c','MarkerFaceColor','c');
plot(ax{9},EV(:,4),EV(:,1),'*m','MarkerFaceColor','m');
plot(ax{10},EV(:,5),EV(:,1),'sk','MarkerFaceColor','k');
%%
%  Find the correct dE
i1 = 1;
for i = 1:length(hb)
    [~,i2] = min(abs(hb(i)-Z));
    for j = i1:i2
        [~,ii(j,1)] = min(abs(dE(i)-dEhor(j,:)));
        dEclose(j,1) = dEhor(j,ii(j));
    end
    i1 = i2+1;
end
plot(ax{4},ii,Z,'*c')
% -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
plot(ax{10},dE,hb,'dk','MarkerFaceColor','k')
plot(ax{11},Ax,hb,'dk','MarkerFaceColor','k')
plot(ax{12},r,hb,'dk','MarkerFaceColor','k')
% -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
ResZmdl = 50;
Zmdl = (ResZmdl:ResZmdl:Z(end))'; Zmdl(end) = Z(end);
AxOut = CLASS_Initials.PossiblePrincipalAxis(PAHV,ao);
plot(ax{4},AxOut.AxMin(:,1),Z,'.k')
plot(ax{4},AxOut.AxMin(:,2),Z,'.k')
%%  -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* GET dE0
MHB = [Z(end)]';
dE0 = CLASS_Initials.PossibledE(Z,Zmdl,AxOut,dEhor);
plot(ax{10},dE0,Zmdl,'sg','MarkerFaceColor','g'); drawnow
%  -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* GET r0
r0 = 0.*ones(size(Zmdl));
plot(ax{12},r0,Zmdl,'sb'); drawnow
% -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* GET Ax0
[Ax0,ax] = FUNC_SpecificLayers(MHB,dE0,r0,Zmdl,ObsDta,dZ,ax);
plot(ax{11},Ax0,Zmdl,'sg','MarkerFaceColor','g'); drawnow
%  -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* TRY IT
ax = FUNC_PlotOptimizedModel(ax,dE0,r0,Ax0,Zmdl,Z,dZ,ao);
%% -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* INITIALIZING INVERRSION
options1 = optimoptions('fmincon','Display','iter-detailed');
options1.UseParallel = true;
options1.OptimalityTolerance = 1e-3;
options1.FunctionTolerance = 1e-3;
options1.StepTolerance = 1;
% -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* INVERSION
resInv = ResZmdl;
Zinv = (resInv:resInv:Z(end))';
% -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* INVERT Ax
% [dE0,r0,Ax,ax] = CLASS_ClassicPointInv.Invert_Ax(dE0,r0,Ax0,Zmdl,Zinv,Z,dZ,ao,ObsDta,options1,ax,[16]);
[dE0,r0,Ax_est,ax] = CLASS_Legendre.Invert_Legendre_Ax(20,dE0,r0,Ax0,Zmdl,Zinv,Z,dZ,ao,ObsDta,options1,ax,[16]);
% -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* INVERT r
[dE0,r_est,Ax_est,ax] = CLASS_ClassicPointInv.Invert_r(dE0,r0,Ax_est,Zmdl,Zinv,Z,dZ,ao,ObsDta,options1,ax,[5]);
% [dE0,r,Ax,ax] = CLASS_Legendre.Invert_Legendre_r(20,dE0,r0,Ax,Zmdl,Zinv,Z,dZ,ao,ObsDta,options1,ax,[5]);
%
[EigVal,dE_end,r_end] = FUNC_GetEigenValues(ax,dE0,r_est,Zmdl,ResZmdl,[]);




















