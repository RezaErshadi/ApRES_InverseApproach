function [EigVal,HA,r] = FUNC_GetEigenValues(ax,HA0,r1,Zinv,dZ,EV)
%     cla(ax{9})
%     plot(ax{9},EV(:,2),EV(:,1),'sy','MarkerFaceColor','y');
%     plot(ax{9},EV(:,3),EV(:,1),'sc','MarkerFaceColor','c');
%     plot(ax{9},EV(:,4),EV(:,1),'sm','MarkerFaceColor','m');
    HA(1,1) = HA0(1);
    plot(ax{10},HA(1),Zinv(1),'.r','LineWidth',2,'MarkerSize',14)
    EigVal(1,:) = EV1stLayer(HA0(1));
    plot(ax{9},EigVal(1,1),Zinv(1),'or','MarkerFaceColor','r')
    plot(ax{9},EigVal(1,2),Zinv(1),'ok','MarkerFaceColor','k')
    plot(ax{9},EigVal(1,3),Zinv(1),'ob','MarkerFaceColor','b')   
    drawnow
    for i = 1:length(Zinv)-1
        % fixed values
        l1T = EigVal(i,1);
        HAT = HA0(i);
        rT = r1(i);
        HAB = HA0(i+1);
        % possible change in case of no solution
        range_r = 1; %3 DC
        range_HA = 0.025; %0.05 DC
        dRange_r = 0.1;
        dRange_HA = 0.005;
        HAlim = [HAB-range_HA:dRange_HA:HAB+range_HA];
        rlim = [rT-range_r:dRange_r:rT+range_r];
        while true
            HAlim(HAlim<0) = [];
            HAlim(HAlim>0.5) = [];
            [~,ix] = min(abs(HAB - HAlim));
            [~,iy] = min(abs(rT - rlim));
            clear('E1')
            clear('g2')
            clear('dist')
            EV1 = EigVal(:,1);
            EV2 = EigVal(:,2);
            EV3 = EigVal(:,3);
%             tic
            parfor j = 1:length(rlim)
                a = rlim(j);
                [t1,t2,t3] = func_getGr2Dist(a,HAlim,l1T,HAT,EV3,dZ,i,j,ix,iy);
                g2(j,:) = t1;
                dist(j,:) = t2;
                E1(j,:) = t3;
            end
%             toc
            nrm_g2 = g2;
            nrm_dist = rescale(dist,0,100);
            %------ best solution
            Aall = nrm_g2+nrm_dist;
            if sum(~isnan(Aall(:)))>0
                break
            else
                rlim = [rlim(1)-1:dRange_r:rlim(1) , rlim(end):dRange_r:rlim(end)+1];
                HAlim = [HAlim(1)-0.025:dRange_HA:HAlim(1) , HAlim(end):dRange_HA:HAlim(end)+0.025];
            end
        end
        [ibst,jbst] = find(Aall==nanmin(Aall(:)));
        ibst = ibst(1);
        jbst = jbst(1);
        %------ Get dE
        HA(i+1,1) = HAlim(jbst);
        plot(ax{10},HA(i+1),Zinv(i+1),'.r','LineWidth',2,'MarkerSize',14)
        %------ Get r
        r(i,1) = rlim(ibst);
        plot(ax{12},r(i),Zinv(i),'.r','LineWidth',2,'MarkerSize',14)  
        %------ Get EigenValues
        EigVal(i+1,1) = E1(ibst,jbst);
        EigVal(i+1,2) = EigVal(i+1,1) + HA(i+1,1);
        EigVal(i+1,3) = 1 - EigVal(i+1,2) - EigVal(i+1,1);
        plot(ax{9},EigVal(i+1,1),Zinv(i+1),'or','MarkerFaceColor','r')
        plot(ax{9},EigVal(i+1,2),Zinv(i+1),'ok','MarkerFaceColor','k')
        plot(ax{9},EigVal(i+1,3),Zinv(i+1),'ob','MarkerFaceColor','b') 
        drawnow
    end
    r(end+1) = r1(end);
    plot(ax{12},r(end),Zinv(end),'.r','LineWidth',2,'MarkerSize',14) 
    drawnow
end

function [g2,dist,E1] = func_getGr2Dist(a,dElim,E1T,dET,EV3,dZ,i,j,ix,iy)
    parfor k = 1:length(dElim)
        b = dElim(k);
        [EV,ErrorF] = dEr2EV(E1T,dET,b,a);
        if ErrorF~= 0
            E1(1,k) = nan;
            LS(1,k) = nan;
            g2(1,k) = nan;
            dist(1,k) = nan;
        else
            E1(1,k) = EV(1);
            LS(1,k) = (EV3(end)-EV(3))/dZ;
            GR1 = diff([EV3 ; EV(3)])./dZ;
            GR2 = diff(GR1)./dZ;
            if i > 1
                g2(1,k) = norm(GR2);
            else
                g2(1,k) = norm(GR1);
            end
            dist(1,k) = sqrt( (ix-k).^2 + (iy-j).^2 );
        end
    end
    LS(abs(LS)>1.5e-3) = nan;
    LS(abs(LS)<1e-6) = nan;
    E1(isnan(LS)) = nan; 
    g2(isnan(LS)) = nan;
    dist(isnan(LS)) = nan;
end

function EigVal1 = EV1stLayer(dE1)
    E1 = 0.333; 
    k = 1;
    while true
        E2 = dE1 + E1;
        E3 = 1-E2-E1;
        ch(1) = E1 >= 0;
        ch(2) = E2>= 0;
        ch(3) = E3>= 0;
        ch(4) = E1 < 0.334;
        ch(5) = E2 <= 0.5;
        ch(6) = E3 > 0.332;
        ch(7) = E3 <= 1;
        ch(8) = E2 >= E1;
        ch(9) = E3 >= E2;
        if sum(ch) == length(ch)
            break
        elseif E1 <= 0
            fprintf("Error, Couldn't find a solution \n")
            break
        else
            E1 = E1-0.00001;
        end
        k = k+1;
    end
    EigVal1(1,:) = sort([E1,E2,E3]);
end
function [EigVal,ErrorF] = dEr2EV(E1T,dET,dEB,rT)
    rT = 10.^(rT./20); % conver from dB to ratio    
    E2T = dET + E1T;
    %----------------- solve the r equation
    syms unknown
    eqn = rT == (((dEB + unknown) - E2T).^2) ./ ((unknown - E1T).^2);
    E1_2 = double(vpasolve(eqn,unknown));    
    %---------------- Check the roots and conditions
    if isempty(E1_2) % the equation has no root
        EigVal = [nan nan nan];
        ErrorF = nan;
    else % the equation has one or more roots
        for i  = 1:length(E1_2)
            EV(i,1) = E1_2(i);
            EV(i,2) = dEB + EV(i,1);
            EV(i,3) = 1-EV(i,2)-EV(i,1);
            ch(i,1) = EV(i,1) >= 0;
            ch(i,2) = EV(i,2)>= 0;
            ch(i,3) = EV(i,3)>= 0;
            ch(i,4) = EV(i,1) < 0.34;
            ch(i,5) = EV(i,2) <= 0.5;
            ch(i,6) = EV(i,3) > 0.32;
            ch(i,7) = EV(i,3)<=1;
            ch(i,8) = EV(i,2) >= EV(i,1);
            ch(i,9) = EV(i,3) >= EV(i,2);
            sch(i) = sum(ch(i,:));
        end
        [~,im] = max(sch);
        ErrorF = length(ch)-sch(im(1));
        if ErrorF == 0
            EigVal = EV(im(1),:);
        else
            EigVal = [nan nan nan];
        end
    end
end