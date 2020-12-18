function [v1_0,ax] = FUNC_SpecificLayers(MHB,HA0,r0,Zmdl,ObsDta,dZ,ax,CostFunc)
%     CP_PM = ObsDta{14};
%     CP_PM(CP_PM>0) = 1;
%     CP_PM(CP_PM<0) = -1;
%     ObsDta{14} = CP_PM;
    v1_0 = [];
    ao1 = 0:5:179;
    ErrorValue = nan(length(MHB),length(ao1));
    s2nr = 0;
    for i = 1:length(MHB)
        [~,i2] = min(abs(MHB(i)-Zmdl));
        cz = Zmdl(1:i2);
        cHA0 = HA0(1:i2);
        cr0 = r0(1:i2);
        AxComb = ao1;
        Ax0Comb = repmat(AxComb,length(cz),1);
    %     ErrorValue = nan(1,size(AxComb,2));
        if i ~= 1
            Ax0Comb(1:length(v1_0),:) = repmat(v1_0,1,length(ao1));
        end 
        ll = length(0:dZ:cz(end));
        for jj = 1:length(ObsDta)
            ObsC{jj} = ObsDta{jj}(1:ll,:);
        end 
        pause(1);
        parfor j = 1:size(AxComb,2)
            OP0 = CLASS_FM.BeginForwardModel(cz,cHA0,cr0,Ax0Comb(:,j),dZ,s2nr);
            EstPar = OP0.Dta;
%             CP_PM2 = EstPar{14};
%             CP_PM2(CP_PM2>0) = 1;
%             CP_PM2(CP_PM2<0) = -1;
%             EstPar{14} = CP_PM2;
%             mf7_temp = FUNC_GetTheMisfit(ObsC,EstPar,7,1); % HH misfit
%             mf14_temp = FUNC_GetTheMisfit(ObsC,EstPar,14,1); % HH misfit
%             mf7(1,j) = norm(mf7_temp);
%             mf14(1,j) = norm(mf14_temp);
            misfit = [];
            for k = 1:length(CostFunc)
                mf = FUNC_GetTheMisfit(ObsC,EstPar,CostFunc(k),1); % HH misfit
                misfit(k) = norm(mf);
            end
            ErrorValue(i,j) = sum(misfit);
        end
        [~,im1] = min(ErrorValue(i,:));
%         cla(ax{11})
        v1_0 = Ax0Comb(:,im1);
%         cla(ax{11})
%         plot(ax{11},Ax0,cz,'sg','MarkerFaceColor','g'); drawnow
    end
end