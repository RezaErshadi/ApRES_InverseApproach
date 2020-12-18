function [v1_0,ax] = FUNC_SpecificLayers3(MHB,HA0,r0,v1_500,Zmdl,ObsDta,dZ,ax,CostFunc)
    v1_0 = v1_500;
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
        Ax0Comb(1:length(v1_0),:) = repmat(v1_0,1,length(ao1));
        ll = length(0:dZ:cz(end));
        for jj = 1:length(ObsDta)
            ObsC{jj} = ObsDta{jj}(1:ll,:);
        end 
        pause(1);
        parfor j = 1:size(AxComb,2)
            OP0 = CLASS_FM.BeginForwardModel(cz,cHA0,cr0,Ax0Comb(:,j),dZ,s2nr);
            EstPar = OP0.Dta;
            misfit = [];
            for k = 1:length(CostFunc)
                mf = FUNC_GetTheMisfit(ObsC,EstPar,CostFunc(k),1); % HH misfit
                misfit(k) = norm(mf);
            end
            ErrorValue(i,j) = sum(misfit);
        end
        [~,im1] = min(ErrorValue(i,:));
        v1_0 = Ax0Comb(:,im1);
    end
end