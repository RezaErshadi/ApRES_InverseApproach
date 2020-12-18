function [Ax0,r0,ax] = FUNC_SpecificLayers2(MHB,HA0,Zmdl,ObsDta,dZ,ax,CostFunc)
    Ax0 = [];
    r0 = [];
    ao1 = 0:5:179;
    r1 = -5 : 1 : 5;
    ErrorValue = nan(length(MHB),length(ao1));
    combinp = combvec(ao1,r1);
    s2nr = 0;
    for i = 1:length(MHB)
        [~,i2] = min(abs(MHB(i)-Zmdl));
        cz = Zmdl(1:i2);
        cHA0 = HA0(1:i2);
        Ax0Comb = repmat(combinp(1,:),length(cz),1);
        r0Comb = repmat(combinp(2,:),length(cz),1);
    %     ErrorValue = nan(1,size(AxComb,2));
        if i ~= 1
            Ax0Comb(1:length(Ax0),:) = repmat(Ax0,1,length(combinp));
            r0Comb(1:length(r0),:) = repmat(r0,1,length(combinp));
        end 
        ll = length(0:dZ:cz(end));
        for jj = 1:length(ObsDta)
            ObsC{jj} = ObsDta{jj}(1:ll,:);
        end 
        pause(1);
        parfor j = 1:size(combinp,2)
            OP0 = CLASS_FM.BeginForwardModel(cz,cHA0,r0Comb(:,j),Ax0Comb(:,j),dZ,s2nr);
            EstPar = OP0.Dta;
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
        Ax0 = Ax0Comb(:,im1);
        r0 = r0Comb(:,im1);
%         cla(ax{11})
%         plot(ax{11},Ax0,cz,'sg','MarkerFaceColor','g'); drawnow
    end
end