% All the related functions to invert HA,r,Ax pointwise
classdef CLASS_ClassicPointInv
    methods(Static)
%%    **********************************************************************
function [HA0,r0,V1,ax] = Invert_v1(HA0,r0,v1_0,Zmdl,Zinv,Z,dZ,ao,ObsDta,options1,ax,CF)
InvPar = ["","Ax"];
all0 = [HA0 r0 v1_0];
i1 = 1;
for i = 1:size(Zinv,1)
    [~,i2] = min(abs(Zinv(i,1)-Zmdl));
    OPI_Ax(i,1) = mean(v1_0(i1:i2));
    i1 = i2+1;
end
bnd_Ax(:,1) = 0.*ones(size(OPI_Ax));
bnd_Ax(:,2) = 180.*ones(size(OPI_Ax));
mp0 = [OPI_Ax];
BND = [bnd_Ax];
tic
fun = @(OPI)CLASS_ClassicPointInv.funfmincon(OPI,all0,Zmdl,Zinv,dZ,ObsDta,InvPar,CF);
OPO = fmincon(fun,mp0,[],[],[],[],BND(:,1),BND(:,2),[],options1);
toc
[HA0,r0,V1] = CLASS_ClassicPointInv.funcGetAll0(OPO,all0,Zmdl,Zinv,InvPar);
V1(V1>179) = V1(V1>179)-180;
V1(V1<0) = v1_0(V1<0)+180;
ax = FUNC_PlotOptimizedModel(ax,HA0,r0,V1,Zmdl,Z,dZ,ao);
plot(ax{11},V1,Zmdl,'.c','MarkerSize',14),drawnow
end
%%
function [HA0,r00,v1_0,ax] = Invert_r(HA0,r0,v1_0,Zmdl,Zinv,Z,dZ,ao,ObsDta,options1,ax,CF)
InvPar = ["r",""];
all0 = [HA0 r0 v1_0];
i1 = 1;
for i = 1:size(Zinv,1)
    [~,i2] = min(abs(Zinv(i,1)-Zmdl));
    OPI_r(i,1) = mean(r0(i1:i2));
    i1 = i2+1;
end
bnd_r(:,1) = -30.*ones(size(Zinv));
bnd_r(:,2) = 30.*ones(size(Zinv));
mp0 = [OPI_r];
BND = [bnd_r];
tic
fun = @(OPI)CLASS_ClassicPointInv.funfmincon(OPI,all0,Zmdl,Zinv,dZ,ObsDta,InvPar,CF);
OPO = fmincon(fun,mp0,[],[],[],[],BND(:,1),BND(:,2),[],options1);
toc
[HA0,r00,v1_0] = CLASS_ClassicPointInv.funcGetAll0(OPO,all0,Zmdl,Zinv,InvPar);
ax = FUNC_PlotOptimizedModel(ax,HA0,r00,v1_0,Zmdl,Z,dZ,ao);
plot(ax{12},r00,Zmdl,'sg','MarkerFaceColor','g'); drawnow
end
%%
function ErrorValue = funfmincon(OPI,all0,Zmdl,Zinv,dZ,ObsDta,InvPar,CF)
    [HA0,r0,v1_0] = CLASS_ClassicPointInv.funcGetAll0(OPI,all0,Zmdl,Zinv,InvPar);
    OP = CLASS_FM.BeginForwardModel(Zmdl,HA0,r0,v1_0,dZ,0);
    EstPar = OP.Dta;
    for i = 1:length(CF)
        msfTyp = 1;
        mf = FUNC_GetTheMisfit(ObsDta,EstPar,CF(i),msfTyp); % HH misfit
        misfit(i) = norm(mf);
    end
    ErrorValue = sum(misfit);
    fprintf("****************************************************** \n");
    for i = 1:length(Zinv)
        [~,i2] = min(abs(Zinv(i)-Zmdl));
        fprintf("HApth: %.0f -------> HA: %.4f && r: %.2f && Ax: %.2f \n",Zinv(i),HA0(i2),r0(i2),v1_0(i2));
    end
    fprintf("****************************************************** \n");
end
%%
function [HA0,r0,v1_0] = funcGetAll0(OPI,all0,Zmdl,Zinv,InvPar)
    HA0 = all0(:,1);
    r0 = all0(:,2);
    v1_0 = all0(:,3);
    %---------------------------
    if InvPar(1) == "r" && InvPar(2) == "Ax"
        OPI_r = OPI(1:length(Zinv));
        OPI_Ax = OPI(length(Zinv)+1:end);
    end
    if InvPar(1) == "r" && InvPar(2) == ""
        OPI_r = OPI(1:length(Zinv));
        OPI_Ax = v1_0;
    end
    if InvPar(2) == "Ax" && InvPar(1) == ""
        OPI_Ax = OPI(1:length(Zinv));
        OPI_r = r0;
    end
    %---------------------------
    if length(Zinv)==length(Zmdl)
        r0 = OPI_r;
        v1_0 = OPI_Ax;
    elseif length(Zinv)~=length(Zmdl)
        i1 = 1;
        for i = 1:length(Zinv)
            [~,i2] = min(abs(Zinv(i)-Zmdl));
            if InvPar(1) == "r" && InvPar(2) == "Ax"
                r0(i1:i2,1) = OPI_r(i);
                v1_0(i1:i2,1) = OPI_Ax(i);
            end
            if InvPar(1) == "r" && InvPar(2) == ""
                r0(i1:i2,1) = OPI_r(i);
            end
            if InvPar(2) == "Ax" && InvPar(1) == ""
                v1_0(i1:i2,1) = OPI_Ax(i);
            end
            i1 = i2+1;
        end
    end
end
%%
end
end