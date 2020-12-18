function [hh,vh,hv,vv,Z,dZ,m_az,f,Bed,MetaData] = FUNC_ReadApRES(DtaDir,maxRange,BedRange)
if ismac || isunix % path seperature
    ps = '/';
elseif ispc
    ps = '\';
end
%%
p=1; % pad factor (i.e. level of interpolation to use during fft)
winFun=@blackman; % window function handle
frange=[2e8 , 4e8]; %Normally [2e8,4e8]
% Density From field observations
rhos=407.0613; % surface density (kg/m^3)
rhoi=907.7165; % ice density (kg/m^3)
Lrho=39.5512;  % half depth decay (m) (assumes rho=rhoi+(rhos-rhoi)*exp(-(H-z)/Lrho))
nI=1.68;
% bulkAlignRange=[80 100];
% dzPowerNorm=10;
% dzdiffphase=10;
% RelPowerCutOff=999999; 
% diffphaseCutOff=99999;
% DepthCutoff=-99999;
%%
% find the MetaData
MetaDataPath = string(DtaDir)+ps+"MetaData.txt";
if isfile(MetaDataPath)
    fid = fopen(MetaDataPath,'r');
    MD1 = fscanf(fid,'%s');
    fclose(fid);
    MD2 = split(MD1,';');
    i = 1;
    while MD2{i} ~= "end"
        a1 = split(MD2{i},':');
        if contains(a1{1},'(')
            a2 = split(a1{1},'(');
            a2 = a2{1};
        else
            a2 = a1{1};
        end
        MetaData(i,1) = string(a2);
        MetaData(i,2) = string(a1{2});
        i = i+1;
    end
else
    MetaData = [];
end
%% find the ApRES files (or fast file)
DataList = [dir(string(DtaDir)+ps+"*.dat") ; dir(string(DtaDir)+ps+"*.DAT") ; dir(string(DtaDir)+ps+"*.mat")];
for i = 1:length(DataList)
    DataInfo(i,1) = string(DataList(i).name);
    DataInfo(i,2) = string(DataList(i).folder);
end
%% Load Data or Saved file
loadtype = 2;
isfast = find(contains(DataInfo(:,1),strcat(string(maxRange),'m_fast.mat')), 1);
if ~isempty(isfast)
        loadtype = 1;
end
switch loadtype
    case 1
        fastpath = strcat(DataInfo(isfast,2),ps,DataInfo(isfast,1));
        load(fastpath);
    case 2
    % find the azimuth and antenna orientation
        DataInfo(~contains(DataInfo(:,1),'deg'),:) = [];
        SplitDataName = split(DataInfo(:,1),'_');
        az = SplitDataName(contains(SplitDataName,'deg'));
        DataInfo(:,3)  = replace(az,'deg','');
        az = str2double(DataInfo(:,3));
        [az,iaz] = sort(az);
        m_az = unique(az'); % measured azimuth
        DataInfo = DataInfo(iaz,:);
        DataInfo(contains(DataInfo(:,1),'HH'),4)  = "HH";
        DataInfo(contains(DataInfo(:,1),'VH'),4)  = "VH";
        DataInfo(contains(DataInfo(:,1),'HV'),4)  = "HV";
        DataInfo(contains(DataInfo(:,1),'VV'),4)  = "VV";
        m_or = unique(DataInfo(:,4)); % measured orientation
        Bed = [];
        hh = []; vh = []; hv = []; vv = [];
        Z = [];
        indNoDta = [];
        f = nan;
        % Load the data
        for i = 1:length(m_az) % read existing azimuths [°]
            for j = 1:length(m_or) % read existing orientations (HH,VH,HV,VV)
                % find all the data (same orientation and azimuth different time)
                iload = find(DataInfo(:,3) == string(m_az(i)) & DataInfo(:,4) == m_or(j));
                [range,specCor,fNew] = func_Concatenate(iload,DataInfo,p,maxRange,winFun,ps,frange);
                if isempty(Z)
                    Z = range;
                    Z = func_TrueDepthCorrection(nI,rhos,rhoi,Lrho,Z); % Correct Depth
                end
                if isnan(f) || fNew~=f
                    f = nanmean([f,fNew]);
                end
                if isnan(specCor)
                    indNoDta(end+1,1) = i;
                else
                    Bed(j,i) = func_Bed(range,specCor,BedRange);
                    Bed(j,i) = func_TrueDepthCorrection(nI,rhos,rhoi,Lrho,Bed(j,i));
                    switch m_or(j)
                        case "HH"
                            hh(:,i) = specCor.';
                        case "VH"
                            vh(:,i) = specCor.';
                        case "HV"
                            hv(:,i) = specCor.';
                        case "VV"
                            vv(:,i) = specCor.';
                    end
                end
            end
        end
        if length(Z) == size(hh,1)
            hh(:,indNoDta) = nan;
        end
        if length(Z) == size(vh,1)
            vh(:,indNoDta) = nan;
        end
        if length(Z) == size(hv,1)
            hv(:,indNoDta) = nan;
        end
        if length(Z) == size(vv,1)
            vv(:,indNoDta) = nan;
        end
        Z(1) = 1e-20; Z = Z'; % depth vector [m]
        dZ = mean(diff(Z)); % depth resolution [m]
        Bed= nanmean(Bed(:));
        %%
        ws_save = split(DtaDir,'/');
        ws_save = strcat(DtaDir,ps,ws_save{end},'_',string(maxRange),'m_fast.mat');
        save(ws_save,'hh','vh','hv','vv','f','Z','dZ','m_az','Bed');
end
end

function [Bed] = func_Bed(range,specCor,BedRange)
    Bed = nan;
    if ~isempty(BedRange)
        bn = fmcw_findbed(range,specCor,BedRange,'xcor');
        Bed = range(bn);
    end
end

function [TrueDepth] = func_TrueDepthCorrection(nI,rhos,rhoi,Lrho,inp)
    TrueDepth = [];
    if ~isempty(inp)
        ShouldBeZero=@(d,dI,L,RhoSp) -dI+d+L*(nI-1)/nI*(1-RhoSp)*(exp(-d/L)-1);
        TrueDepthfun=@(RhoSp,L,dI) FUNC_bisection(@(d) ShouldBeZero(d,dI,L,RhoSp),0,max(dI)+20);
        TrueDepth = TrueDepthfun(rhos/rhoi,Lrho,inp);
    end
end

function [specCor,ax] = func_Signal2Noise(specCor,range,maxRange,BedDepth,trshld,plt)
    SignaldB = abs(specCor);
    NoiseFloor=max(SignaldB(range>maxRange-trshld));
    NoiseFloordb=20*log10(NoiseFloor);
    ax = [];
    if plt == 1
        figure; hold all;
        ax = gca;
        plot(20*log10(abs(specCor)),range,'.r');
        set(ax,'Ydir','reverse')
        set(ax,'Color',[0.6 0.6 0.6]);
        if ~isempty(BedDepth)
            plot([NoiseFloordb,NoiseFloordb],[0,BedDepth],'-b','linewidth',2);
            plot([NoiseFloordb,NoiseFloordb],[BedDepth maxRange],'--b','linewidth',2);
            plot([min(20*log10(SignaldB)),max(20*log10(SignaldB))],[BedDepth,BedDepth],'--k','linewidth',2);
        else
            plot([NoiseFloordb,NoiseFloordb],[0,maxRange],'--k','linewidth',2);
        end
        plot(20*log10(abs(specCor(abs(specCor)>NoiseFloor))),range(abs(specCor)>NoiseFloor),'.g');
        title('Noise | Signal')
        xlabel('Power [dB]')
        ylabel('Depth [m]')
        grid on
        xlim([min(20*log10(SignaldB)),max(20*log10(SignaldB))])
        ylim([range(1) range(end)])
    end
    specCor(SignaldB<NoiseFloor)=nan;
end

function [range,specCor,f] = func_Concatenate(iload,DataInfo,p,maxRange,winFun,ps,frange)
    range0 = [];
    specCor0 = [];
    f0 = [];
    for k=1:length(iload) % average in case of several measurements with same setup (time laps)
        filePath = strcat(DataInfo(iload(k),2),ps,DataInfo(iload(k),1));
        try
            DtaLoad = fmcw_load(filePath,1); % load the ApRES file
            DtaLoad=fmcw_cull_freq(DtaLoad,frange);
            if size(DtaLoad.vif,1) > 1
                % Takes mean of fmcw_burst
                DtaMean = fmcw_burst_mean(DtaLoad); 
            else
                DtaMean = DtaLoad;
            end
            [TempRange,~,TempSpecCor,~] = fmcw_range(DtaMean,p,maxRange,winFun);
            range0(end+1,:) = TempRange;
            specCor0(end+1,:) = TempSpecCor; 
            f0(end+1) = mean(DtaMean.f); % center frequency
        catch
            fprintf("Error in loading file: " + DataInfo(iload(k),1) +"\n")
        end
    end
    if ~isempty(range0)
        range = mean(range0,1);
        specCor = mean(specCor0,1);
        f = nanmean(f0);
    else
        range = nan;
        specCor = nan;
        f = nan;
    end
end





