classdef CLASS_Denoising
    methods(Static)
%% Run Denoising
function [PA,PD] = RunDenoising(DenoisingFlag,PA,PD,dZ)
    if ~isempty(DenoisingFlag)
        [PA] = CLASS_Denoising.DenoiseParam(PA,dZ,DenoisingFlag.PA);
        [PD] = CLASS_Denoising.DenoiseParam(PD,dZ,DenoisingFlag.PD);
    end
end
%% Moving Average
    function [otp] = MovingAverage(inp,NoiseWin,dZ)
        stwnd = NoiseWin./dZ;
        otp = inp;
        otp = smoothdata(otp,1,'movmean',[stwnd],'omitnan');
    end
%% Average Depth
     function otp=AverageDepth(inp,Z,Zwin)
         otp = nan(size(inp));
        for i=1:length(Z)
            iok=Z>Z(i)-Zwin/2&Z<Z(i)+Zwin/2;
            otp(i,:)=mean(inp(iok,:),1);
        end
     end   
%% 1D Convolution
    function otp = Conv1D(inp,ConvWin,dZ)
            win=gausswin(round(ConvWin/dZ,0));
            win=win/sum(win);
            for i = 1:size(inp,2)
                otp(:,i) = conv(inp(:,i),win,'same');
            end
    end
%% 2D Convolution
    function otp = Conv2D(inp,ConvWin,dZ)
            win=gausswin(round(ConvWin/dZ,0));
            win=win/sum(win);
            otp = conv2(inp,win,'same');
    end
%% PCA
    function otp = DenoisePCA(inp,slct_eig)
        %% demean
        inp(isinf(inp)) = nan;
        inp = fillmissing(inp,'nearest');
        nv = size(inp,2); % number of vectores (data)
        X = inp;
        avg = nanmean(X,1);
        B = X - avg;
        %%
        c = B' * B; % covariance of the demean
        %%
        [V,D,W] = eig(c);
        % D = diagonal matrix of eigen values
        % W = left eifen vectors: W'*c = D*W'
        % V = right eifen vectors: c*V = V*D
        [~,I] = sort(diag(D),'descend');
        Vs = V(:,I); % c*Vs = V*D
        %%
        T = B * Vs;
        for i=1:nv
            nrmT(:,i)=T(:,i)/norm(T(:,i),2); % normalizing EIGEN_DIGITS
        end
        %%
        otp=zeros(size(B));
        sEV =slct_eig;
        for k=1:nv
            unq=[];
            for i=1:length(sEV)
                tmp = (nrmT(:,sEV(i)))'*(B(:,k))*(nrmT(:,sEV(i)));
                unq=[unq  tmp];
            end
            otp(:,k)=sum(unq,2)+avg(k);
        end
    end
%% DenoiseParameter
    function [inp] = DenoiseParam(inp,dZ,DenFlag)
        DenFlag(DenFlag(:,1)=="0",:) = [];
        [~,isrt] = sort(str2double(DenFlag(:,1)));
        DenFlag = DenFlag(isrt,:);
        for i = 1:size(DenFlag,1)
            switch DenFlag(i,2)
                case "MovingAverage"
                    NoiseWin = str2double(DenFlag(i,3));
                    [inp.HH] = CLASS_Denoising.MovingAverage(inp.HH,NoiseWin,dZ);
                    [inp.VH] = CLASS_Denoising.MovingAverage(inp.VH,NoiseWin,dZ);
                    [inp.HV] = CLASS_Denoising.MovingAverage(inp.HV,NoiseWin,dZ);
                    [inp.VV] = CLASS_Denoising.MovingAverage(inp.VV,NoiseWin,dZ);
                case "Conv1D"
                    ConvWin = str2double(DenFlag(i,3));
                    [inp.HH] = CLASS_Denoising.Conv1D(inp.HH,ConvWin,dZ);
                    [inp.VH] = CLASS_Denoising.Conv1D(inp.VH,ConvWin,dZ);
                    [inp.HV] = CLASS_Denoising.Conv1D(inp.HV,ConvWin,dZ);
                    [inp.VV] = CLASS_Denoising.Conv1D(inp.VV,ConvWin,dZ);
                case "Conv2D"
                    ConvWin = str2double(DenFlag(i,3));
                    [inp.HH] = CLASS_Denoising.Conv1D(inp.HH,ConvWin,dZ);
                    [inp.VH] = CLASS_Denoising.Conv1D(inp.VH,ConvWin,dZ);
                    [inp.HV] = CLASS_Denoising.Conv1D(inp.HV,ConvWin,dZ);
                    [inp.VV] = CLASS_Denoising.Conv1D(inp.VV,ConvWin,dZ);  
                case "DenoisePCA"
                    slct_eig = str2double(DenFlag(i,3));
                    [inp.HH] = CLASS_Denoising.DenoisePCA(inp.HH,slct_eig);
                    [inp.VH] = CLASS_Denoising.DenoisePCA(inp.VH,slct_eig);
                    [inp.HV] = CLASS_Denoising.DenoisePCA(inp.HV,slct_eig);
                    [inp.VV] = CLASS_Denoising.DenoisePCA(inp.VV,slct_eig);
            end
        end  
    end
    end
end