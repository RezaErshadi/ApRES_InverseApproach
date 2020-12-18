classdef CLASS_InvPlot
    methods(Static)
%%  
    function [fg,ax] = InversionPlot(HH,HV,CP,Psi,Z,ao,bd,figvis)
        %%
        pltdim = [0.5000    0.0256    0.5003    0.9213];
        sbplt_row = 3;
        sbplt_col  = 4;
        CLASS_FixedPlot.SetFigureSize(pltdim(1),pltdim(2),pltdim(3),pltdim(4));    fg = gcf;   fg.Color = "white";
        if figvis ~= 1
            set(fg,'Visible','off');
        end
        cm1 = getPyPlot_cMap('seismic',100);
%         cm1 = CLASS_FixedPlot.redblue(100);
%         cm2 = copper(100);
%         cm21 = jet(100);
%         cm2(101:200,:) = cm21;
%         cm2(end,:) = [1 1 1];
%         cm2(1,:) = [1 1 1];
%         cm2(99,:) = [1 1 1];
%         cm2(100,:) = [1 1 1];
%         cm2(101,:) = [1 1 1];
%         cm2(102,:) = [1 1 1];
%         cm3 = CLASS_FixedPlot.redblue(100);
%         cm3(1,:) = [0 0 0];
%         cm3(end,:) = [0 0 0];
        ll = -0.01; % girdle strength lower limit
        cmlength = length(-0.5:abs(ll):0.5);
        Psi(Psi<ll) = nan;
        Psi(Psi>0.5) = nan;    
        cm3 = getPyPlot_cMap('tab20b',cmlength-1);
%         cm3 = summer(cmlength-1);
        fnt = 14;
        limY = [Z(1) Z(end)];
        yl = linspace(0,Z(end)-mod(Z(end),5),5);
        yl(1) = Z(1);
        cgcl = [0.9 0.9 0.9];
        xl_theta = 0:45:180;
        %%
        %------------------------------------
        iax = 1;
        ax{iax} = subplot(sbplt_row,sbplt_col,1);
        imagesc(ax{iax},ao,Z,HH)
        colormap(ax{iax},cm1)
        cb1 = colorbar(ax{iax},'Location','southoutside');
        caxis(ax{iax},[-5 5])
        set(ax{iax},'YDIR','reverse')
        ylim(ax{iax},limY)
        yticks(ax{iax},yl);
        yticklabels(ax{iax},["0" string(yl(2:end))]);
        xlim(ax{iax},[0 180]);
%         xticks(ax{iax},xl_theta)
        ylabel(ax{iax},'Depth [m]');
%         xlabel(ax{iax},"\alpha [deg]");
        xticklabels(ax{iax},[])
        set(ax{iax},'FontSize',fnt)
        title(ax{iax},['\fontsize{14}(a) HH Power Anomaly'])
        grid(ax{iax},'on')
        hold(ax{iax},'on');
        if isempty(bd) ~= 1
            plot(ax{iax},xlim(ax{iax}),[bd bd],'-k','LineWidth',3);
        end
        %------------------------------------
        iax = 2;
        ax{iax} = subplot(sbplt_row,sbplt_col,2);
        imagesc(ax{iax},ao,Z,CP)
        colormap(ax{iax},cm1)
        cb2 = colorbar(ax{iax},'Location','southoutside');
        caxis(ax{iax},[-pi pi])
        set(ax{iax},'YDIR','reverse')
        ylim(ax{iax},limY)
        yticks(ax{iax},yl);
        yticklabels(ax{iax},[])
        xlim(ax{iax},[0 180]);
%         xticks(ax{iax},xl_theta)
%         xlabel(ax{iax},"\alpha [deg]");
        xticklabels(ax{iax},[])
        set(ax{iax},'FontSize',fnt)
        title(ax{iax},['\fontsize{14}(b) HHVV Coherence Phase'])
        grid(ax{iax},'on')
        hold(ax{iax},'on');
        if isempty(bd) ~= 1
            plot(ax{iax},xlim(ax{iax}),[bd bd],'-k','LineWidth',3);
        end
        %------------------------------------
        iax = 3;
        ax{iax} = subplot(sbplt_row,sbplt_col,3); 
        imagesc(ax{iax},ao,Z,HV)
        colormap(ax{iax},cm1)
        cb3 = colorbar(ax{iax},'Location','southoutside');
        caxis(ax{iax},[-5 5])
        set(ax{iax},'YDIR','reverse')
        ylim(ax{iax},limY)
        yticks(ax{iax},yl);
        yticklabels(ax{iax},[])
        xlim(ax{iax},[0 180]);
%         xticks(ax{iax},xl_theta)
%         xlabel(ax{iax},"\alpha [deg]");
        xticklabels(ax{iax},[])
        set(ax{iax},'FontSize',fnt)
        title(ax{iax},['\fontsize{14}(c) HV Power Anomaly'])
        grid(ax{iax},'on')
        hold(ax{iax},'on');
        if isempty(bd) ~= 1
            plot(ax{iax},xlim(ax{iax}),[bd bd],'-k','LineWidth',3);
        end
        %------------------------------------
        iax = 4;
        ax{iax} = subplot(sbplt_row,sbplt_col,4);
        imAlpha=ones(size(Psi));
        imAlpha(isnan(Psi))=0;
        imagesc(ax{iax},ao,Z,Psi,'AlphaData',imAlpha)
        colormap(ax{iax},cm3)
        cb4 = colorbar(ax{iax},'Location','southoutside');
        caxis(ax{iax},[ll 0.5])
        set(ax{iax},'YDIR','reverse')
        ylim(ax{iax},limY)
        yticks(ax{iax},yl);
        yticklabels(ax{iax},[])
        xlim(ax{iax},[0 180]);
%         xticks(ax{iax},xl_theta)
%         xlabel(ax{iax},"\alpha [deg]");
        xticklabels(ax{iax},[])
        set(ax{iax},'FontSize',fnt)
        title(ax{iax},['\fontsize{14}(d) Scaled Phase Derivative'])
        grid(ax{iax},'on')
        hold(ax{iax},'on');
        if isempty(bd) ~= 1
            plot(ax{iax},xlim(ax{iax}),[bd bd],'-k','LineWidth',3);
        end
        %------------------------------------
        iax = 5;
        ax{iax} = subplot(sbplt_row,sbplt_col,5);
        %------------------------------------
        iax = 6;
        ax{iax} = subplot(sbplt_row,sbplt_col,6);
        %------------------------------------
        iax = 7;
        ax{iax} = subplot(sbplt_row,sbplt_col,7);
        %------------------------------------
        iax = 8;
        ax{iax} = subplot(sbplt_row,sbplt_col,8);
        %------------------------------------
        iax = 9;
        ax{iax} = subplot(sbplt_row,sbplt_col,9);
        set(ax{iax},'Color',cgcl);
        set(ax{iax},'YDIR','reverse')
        ylim(ax{iax},limY)
        yticks(ax{iax},yl);
        yticklabels(ax{iax},["0" string(yl(2:end))]);
        xlim(ax{iax},[0 1])
        xlim(ax{iax},[0 1]);
        xticks(ax{iax},[0:0.25:1]);
        ylabel(ax{iax},'Depth [m]');
        xlabel(ax{iax},"\lambda [-]");
        set(ax{iax},'FontSize',fnt)
        title(ax{iax},['\fontsize{14}(i) Eigenvalues'])
        grid(ax{iax},'on')
        hold(ax{iax},'on');
        %------------------------------------
        iax = 10;
        ax{iax} = subplot(sbplt_row,sbplt_col,10);
        set(ax{iax},'Color',cgcl);
        set(ax{iax},'YDIR','reverse')
        ylim(ax{iax},limY)
        yticklabels(ax{iax},[])
        yticks(ax{iax},yl); 
        xlim(ax{iax},[0 0.50001])
        xticks(ax{iax},[0:0.1:0.5]);
        xlabel(ax{iax},"\Delta\lambda [-]"); 
        set(ax{iax},'FontSize',fnt)
        title(ax{iax},['\fontsize{14}(j) Horizontal Anisotropy'])
        grid(ax{iax},'on')
        hold(ax{iax},'on');
        %------------------------------------
        iax = 11;
        ax{iax} = subplot(sbplt_row,sbplt_col,11);
        set(ax{iax},'Color',cgcl);
        set(ax{iax},'YDIR','reverse')
        ylim(ax{iax},[Z(1) Z(end)])
        yticklabels(ax{iax},[])
        yticks(ax{iax},yl);
        xlim(ax{iax},[0 180])
        xticks(ax{iax},xl_theta);
        xlabel(ax{iax},"\alpha [deg]"); 
        set(ax{iax},'FontSize',fnt)
        title(ax{iax},['\fontsize{14}(k) Horizontal Eigenvectors'])
        grid(ax{iax},'on')
        hold(ax{iax},'on');
        %------------------------------------
        iax =12;
        ax{iax} = subplot(sbplt_row,sbplt_col,12); 
        set(ax{iax},'Color',cgcl);
        set(ax{iax},'YDIR','reverse')
        ylim(ax{iax},limY)
        yticks(ax{iax},yl);
        yticklabels(ax{iax},[])
        xlim(ax{iax},[-30 30])
        xticks(ax{iax},[-30:15:30])
        xlabel(ax{iax},"r [dB]");
        set(ax{iax},'FontSize',fnt)
        title(ax{iax},['\fontsize{14}(l) Reflection Ratio'])
        grid(ax{iax},'on')
        hold(ax{iax},'on');
%         set(ax{iax},'xscale','log')
%         xlim(ax{iax},[1e-3 1e3])
%         xticks(ax{iax},[1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3])
        %%
        pcb1 = get(cb1,'position');
        pcb2 = get(cb2,'position');
        pcb3 = get(cb3,'position');
        pcb4 = get(cb4,'position');
        d = 0.045;
        set(cb1,'Position',[pcb1(1) pcb1(2)-d pcb1(3) pcb1(4)]);
        set(cb2,'Position',[pcb2(1) pcb2(2)-d pcb2(3) pcb2(4)]);
        set(cb3,'Position',[pcb3(1) pcb3(2)-d pcb3(3) pcb3(4)]);
        set(cb4,'Position',[pcb4(1) pcb4(2)-d pcb4(3) pcb4(4)]);
        set(fg,'Position',pltdim);
        cb1.Ticks = [-5 0 5];
%         cb1.TickLabels = string(cb1.Ticks);
        cb2.Ticks = [round(-pi,2) 0 round(pi,2)];
%         cb2.TickLabels = [cb2.Ticks];
        cb3.Ticks = [-5 0 5];
%         cb3.TickLabels = string(cb3.Ticks);
        cb4.Ticks = [0 0.1 0.2 0.3 0.4 0.5];
%         cb4.TickLabels = string(cb4.Ticks);
        
        d = 0.1;
        cb1.Label.String = 'P_{HH} [dB]';
        cb1.Label.FontWeight = 'bold';
        cb1.Label.FontSize = 10;
        Lcb1 = get(cb1.Label,'position');
%         cb1.Label.Rotation = 0;
        cb1.Label.Position = [Lcb1(1) Lcb1(2)+d Lcb1(3)];

        cb2.Label.String = '\phi_{HHVV} [rad]';
        cb2.Label.FontWeight = 'bold';
        cb2.Label.FontSize = 10;
        Lcb2 = get(cb1.Label,'position');
%         cb1.Label.Rotation = 0;
        cb2.Label.Position = [Lcb2(1) Lcb2(2)+d Lcb2(3)];
        
        cb3.Label.String = 'P_{HV} [dB]';
        cb3.Label.FontWeight = 'bold';
        cb3.Label.FontSize = 10;
        Lcb3 = get(cb3.Label,'position');
%         cb1.Label.Rotation = 0;
        cb3.Label.Position = [Lcb3(1) Lcb3(2)+d Lcb3(3)];
        
        cb4.Label.String = '\Psi [-]';
        cb4.Label.FontWeight = 'bold';
        cb4.Label.FontSize = 10;
        Lcb4 = get(cb4.Label,'position');
%         cb1.Label.Rotation = 0;
        cb4.Label.Position = [Lcb4(1) Lcb4(2)+d Lcb4(3)];
        
    txt1 = text('String','Radar Observations','Position',[-297.143942484911 -3566.6686251267 0]);
    set(txt1,'FontSize',18)
    set(txt1,'FontWeight','Bold')
    set(txt1,'Rotation',90)

    txt2 = text('String','Optimized Models','Position',[-297.143942484911 -900.592018742136 0]);
    set(txt2,'FontSize',18)
    set(txt2,'FontWeight','Bold')
    set(txt2,'Rotation',90)
    
    txt3 = text('String','Fabric Parameters','Position',[-297.143942484911 1899.74743832366 0]);
    set(txt3,'FontSize',18)
    set(txt3,'FontWeight','Bold')
    set(txt3,'Rotation',90)
        
        
    end
%%
    function [ax] = UpdateInversionPlot(ax,HH,HV,PCP,dE,Z,cZ,ao,MP,EigVal,figvis)
        cm1 = getPyPlot_cMap('seismic',100);
        cm2 = copper(40);
        cm21 = jet(40);
        cm2(41:80,:) = cm21;
        cm2(end,:) = [1 1 1];
        cm2(1,:) = [1 1 1];
%         cm3 = CLASS_FixedPlot.redblue(100);
%         cm3(1,:) = [0 0 0];
%         cm3(end,:) = [0 0 0];
        ll = -0.01; % girdle strength lower limit
        cmlength = length(-0.5:abs(ll):0.5);
        dE(dE<ll) = nan;
        dE(dE>0.5) = nan;    
        cm3 = getPyPlot_cMap('tab20b',cmlength-1);
        fnt = 14;
        limY = [Z(1) Z(end)];
        yl = linspace(0,Z(end)-mod(Z(end),5),5);
        yl(1) = Z(1);
        xl_theta = 0:45:180;
        %------------------------------------
        if ~isempty(EigVal)
            E1 = EigVal(1);
            E2 = EigVal(2);
            E3 = EigVal(3);
        end
        %------------------------------------
        iax = 5;
        imagesc(ax{iax},ao,Z,HH)
        set(ax{iax},'YDIR','reverse')
        caxis(ax{iax},[-5 5])
        colormap(ax{iax},cm1)
        ylim(ax{iax},limY)
        yticks(ax{iax},yl);
        yticklabels(ax{iax},["0" string(yl(2:end))]);
        xlim(ax{iax},[0 180]);
        xticks(ax{iax},xl_theta)
        ylabel(ax{iax},'Depth [m]');
        xlabel(ax{iax},"\alpha [deg]"); 
        set(ax{iax},'FontSize',fnt)
        title(ax{iax},['\fontsize{14}(e) HH Power Anomaly'])
        grid(ax{iax},'on')
        %------------------------------------
        iax = 6;
        imagesc(ax{iax},ao,Z,PCP)
        set(ax{iax},'YDIR','reverse')
        caxis(ax{iax},[-pi pi])
        colormap(ax{iax},cm1)
        ylim(ax{iax},limY)
        yticks(ax{iax},yl);
        yticklabels(ax{iax},[])
        xlim(ax{iax},[0 180]);
        xticks(ax{iax},xl_theta)
        xlabel(ax{iax},"\alpha [deg]"); 
        set(ax{iax},'FontSize',fnt)
        title(ax{iax},['\fontsize{14}(f) HHVV Coherence Phase'])
        grid(ax{iax},'on')
        %------------------------------------
        iax = 7;
        imagesc(ax{iax},ao,Z,HV)
        set(ax{iax},'YDIR','reverse')
        caxis(ax{iax},[-5 5])
        colormap(ax{iax},cm1)
        ylim(ax{iax},limY)
        yticks(ax{iax},yl);
        yticklabels(ax{iax},[])
        xlim(ax{iax},[0 180]);
        xticks(ax{iax},xl_theta)
        xlabel(ax{iax},"\alpha [deg]"); 
        set(ax{iax},'FontSize',fnt)
        title(ax{iax},['\fontsize{14}(g) HV Power Anomaly'])
        grid(ax{iax},'on')
        %------------------------------------
        iax = 8;
        imAlpha=ones(size(dE));
        imAlpha(isnan(dE))=0;  
        imagesc(ax{iax},ao,Z,dE,'AlphaData',imAlpha)
        set(ax{iax},'YDIR','reverse')
        caxis(ax{iax},[ll 0.5])
        colormap(ax{iax},cm3)
        ylim(ax{iax},limY)
        yticklabels(ax{iax},[])
        yticks(ax{iax},yl);
        xlim(ax{iax},[0 180]);
        xticks(ax{iax},xl_theta)
        xlabel(ax{iax},"\alpha [deg]"); 
        set(ax{iax},'FontSize',fnt)
        title(ax{iax},['\fontsize{14}(h) Scaled Phase Derivative'])
        grid(ax{iax},'on')
        %------------------------------------
        iax = 9;
        if ~isempty(EigVal)
            plot(ax{iax},E1,cZ(end),'or','MarkerFaceColor','r')
            plot(ax{iax},E2,cZ(end),'ok','MarkerFaceColor','k')
            plot(ax{iax},E3,cZ(end),'ob','MarkerFaceColor','b')
        end
        %------------------------------------
        iax = 10;
        if ~isempty(MP) 
            plot(ax{iax},MP(1),cZ(end),'.r','LineWidth',2,'MarkerSize',12)
        end
        %------------------------------------
        iax = 11;
        if ~isempty(MP) 
            plot(ax{iax},MP(3),cZ(end),'.r','LineWidth',2,'MarkerSize',12)
        end
        %------------------------------------
        iax = 12;
        if ~isempty(MP) 
            plot(ax{iax},MP(2),cZ(end),'.r','LineWidth',2,'MarkerSize',12)
%             semilogx(ax{iax},MP(2),cZ(end),'.r','LineWidth',2,'MarkerSize',12)
        end
        %------------------------------------
        if figvis == 1
            drawnow
        end
    end
%%     
    end
end