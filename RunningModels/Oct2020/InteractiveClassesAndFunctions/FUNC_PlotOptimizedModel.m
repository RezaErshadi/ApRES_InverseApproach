function [ax,OP0] = FUNC_PlotOptimizedModel(ax,HA,r,Ax,Zinv,Z,dZ,ao)
    for i=5:8
        cla(ax{i});
    end
    tic
    OP0 = CLASS_FM.BeginForwardModel(Zinv,HA,r,Ax,dZ,0);
    toc
    ObsDta = OP0.Dta;
    estPAHH = ObsDta{5};
    estPAHV = ObsDta{7};
    estCP = ObsDta{14};
    estPsi = ObsDta{18};
    ax = CLASS_InvPlot.UpdateInversionPlot(ax,estPAHH,estPAHV,estCP,estPsi,Z,[],ao,[],[],1);
end