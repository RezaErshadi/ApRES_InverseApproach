function fmcw_ts(cfgfile)

% fmcw_ts(cfgfile)
%
% Process sequence of FMCW shots to calculate timeseries of melt rates
%
% intput:
% cfgfile = text file of settings and input data filenames
% this includes internal and bed ranges and files to process process.
% If this is absent the user is prompted to make one.
%
% Output: matfile with the same name as the cgffile.
%
% Processing steps
% (1) define a cfg file which lists the files to process, the chosen
% internal and bed depths, vertical strain rate for the site and metadata.
% (2) sequentially load and phase process each burst (sub-function getSpec)
%       (a) load burst
%       (b) pre-process: split by attenuator
%                        keep only selected chirp numbers (from cfg)
%       (c) clean: frequency cull, remove bad chirps, remove noisey chirps
%       (d) phase process for spectrum
% (3) determine bed depth 
% (4) Phase calibration and burst stats (subfunction derampStats)
%       (a) for each chirp spectrum remove phase drift through the burst 
%           by differencing the phase in each bin with the bed phase
%       (b) calculate burst stats (mean and standard error in each bin)
% (5) For each pair of burst means: F and G, calculate the inter-burst phase difference 
%       in each bin as pd = ang(FG), where FG = conj(F).*G
% (6) Estimate the inter-burst phase difference from selected internals by either:
%       (a) vector sum: 
%           PD = ang(sum(FG)) (over the selected internals only), or
%       (b) weighted least squares fit:
%           PD = atan(lscov(real(FG),imag(FG),1./(se.^2)))
%           where se is the standard error within each bin, which reduces
%           the weighting applied to bins with high standard error
% (7) Convert the phase change to range and divide by the inter-burst time
%       interval to calculate a thinning rate
% (8) Using the user defined strain rate calculate the strain contribution 
% to thinning and remove this to calculate meltrate
% (9) Save output to .mat file.
% (10) Call external plotting function

% Craig Stewart
% 2014/11/15

% Constants
daysPerYear = 365.25;

% Select cfg file if not entered
if nargin == 0
    cfgfile = uigetfile('*.cfg','Select a cfg file - or cancel to create new');
end
if cfgfile==0
    cfgfile = fmcw_batchlist_write2; % Create cfg file
end
[cfg,file] = fmcw_batchlist_read2(cfgfile); % Read cfg file

%% Processing options
% Pre-processing

% Attenuator settings to keep
cfg.attSetKeep = 'min'; %  'all','min' % all or minimum attenuation
% Chirp selection (for imaging system)
cfg.chirpNo = []; % chirps to use (empty for all)
if ~isempty(cfg.chirpNo)
    % Warn
    disp(['Using chirps: ' mat2str(cfg.chirpNo) ' only'])
end

% Cleaning
%cfg.prepro.freqUseRange = [200.2e6 400e6];
%cfg.prepro.freqUseRange = [201e6 399e6];
cfg.prepro.freqUseRange = [200e6 400e6];
cfg.prepro.cullFirstChirps = 0;
cfg.prepro.keepCleanChirps = 6;

% Internal selection
cfg.bedSeachMethod = 'maxAmp';
cfg.intMode = 'peaks'; % 'all','peaks','uninterp'
cfg.intMinPeakSep = 2; % minimmum separation between internal peaks

% Phase processing settings
cfg.padfactor = 1;

% Phase difference method
cfg.intCombMethod = 'errorWeightFit'; % 'ampWeightSum' 'errorWeightFit' 'meanphasediff'
cfg.errorMode = 'btwnBin'; % 'inBin','btwnBin','combo' (only used if cfg.intCombMethod = 'ampWeightSum')
cfg.leapN = 1; % leapfrog level - difference with non-adjacent bursts
% i.e. compare burst(n) with burst(n+cfg.leapN)

% Plot options
cfg.doPlotTest = 1; % plotting intermediate data for debugging
cfg.plotCI = [1:2];

disp(' ')
disp(['Processing using config file: ' cfgfile])
disp('settings:')
disp(cfg)
disp(' ')

%% Count bursts so we can preallocate matricies before processing
nbursts = zeros(1,length(file));
disp('Checking for number of bursts to process')
for fn = 1:length(file)
    filename = char(file(fn).name);
    disp(['file ' int2str(fn) '/' int2str(length(file)) ':' filename])
    if isnan(file(fn).burst) % use all bursts from file
        nbursts(fn) = fmcw_nbursts(filename);
    end
end
N = sum(nbursts);
[I,Ie,Ic,Imd,bedInd,timeC,dt,PD,PDe] = deal(ones(N-cfg.leapN,1)); % differences

%% Loop through all files and bursts averaging all chirps within the burst
bi = 0; % burst index
ci = 0; % comparison index (always bi-1)
for fn = 1:length(file)
    filename = char(file(fn).name);
    for burst = 1:nbursts(fn)
        bi = bi+1; % master index of burst
        if bi <= cfg.leapN
            %% Load burst
            R(bi) = getSpec(filename,burst,cfg); % load spectrum into structure of raw data
            % Find bed
            bedInd = fmcw_findbed(R(bi).rangeCoarse,abs(mean(R(bi).specRaw,1)),cfg.bedSearchRange,cfg.bedSeachMethod); % bed ind
            S(bi) = derampStats(R(bi),bedInd); % stats structure
            if bi == 1
                % Find internals (just once - assuming these don't move range
                % bins throughout deploy - this may need revision for very long
                % deploys if using "peaks")
                rangeCoarse = R(1).rangeCoarse; % lets assume this doesn't change throughout the timeseries

                switch cfg.intMode
                    case 'all'
                        isInt = rangeCoarse>=min(cfg.intRange) & rangeCoarse<=max(cfg.intRange); % is in internal range
                    case 'uninterp'
                        isInt = rangeCoarse>=min(cfg.intRange) & rangeCoarse<=max(cfg.intRange)...
                            & ~mod([1:length(rangeCoarse)],cfg.padfactor); % is in internal range and is not an interpolated point
                    case 'peaks'
                        % Now only select peaks to get only best SNR
                        dr = mean(diff(rangeCoarse));
                        mpdb = round(cfg.intMinPeakSep/dr); % min peak separation distance (in bins)
                        [~,pkind] = findpeaks(abs(S(1).mean),'MinPeakDistance',mpdb); % amplitude peaks separated by 4m
                        isPeak = zeros(size(S(1).mean));
                        isPeak(pkind) = 1;
                        isInt = rangeCoarse>=min(cfg.intRange) & rangeCoarse<=max(cfg.intRange) & isPeak; % is in internal range
                end
                %intDepth = rangeCoarse(isInt);
                
                % Preallocate large matricies
                %[FG,FGSEP] = deal(ones(N-1,numel(rangeCoarse))); % differences
                %[FGI] = deal(ones(N-1,sum(isInt))); % differences internals only
            end
            
        else
            %% Compare shots
            ci = ci + 1; % comparison index;
            Fi = 1 + mod(ci-1,cfg.leapN); % buffer index
            
            % Check whether this is the bed it was loaded with (they can differ as the bed depth is 
            % set to match the first burst in the comparison
            bedInd(ci) = fmcw_findbed(rangeCoarse,abs(S(Fi).mean),cfg.bedSearchRange,cfg.bedSeachMethod); % bed ind
            if bedInd(ci) == S(Fi).bedInd
                F = S(Fi); % Get burst F directly from burst stats buffer (correct bed applied when loaded)
            else
                % F was previously processed with a non-local bed, so reprocess F from raw data using its own bed location
                F = derampStats(R(Fi),bedInd(ci)); % calculate stats using bedInd of F
            end
                        
            % Load burst G into buffer overwriting F that we've just used
            R(Fi) = getSpec(filename,burst,cfg); % load new burst
            S(Fi) = derampStats(R(Fi),bedInd(ci)); % calculate stats using bedInd of F
            G = S(Fi);
            
            % Comparison properties
            timeC(ci) = mean([F.time G.time]);
            dt(ci) = G.time-F.time;
            disp([int2str(ci) ') F: ' F.filename ' b:' int2str(F.burst) ' t:' datestr(F.time) ' and G: ' G.filename ' b:' int2str(G.burst) ' t:' datestr(G.time)])
            
            %% Compare bursts
            fg = conj(F.mean).*G.mean;
            fgsep = sqrt(F.stePhase.^2 + G.stePhase.^2); % standard error (phase) of fg (assuming uncorrelated noise)
            fgse = sqrt(2)*fgsep.*abs(fg); % standard error of fg
            I(ci) = sum(fg(isInt)); % sum fg over internals
            Ic(ci) = I(ci)./(sqrt(sum(abs(F.mean(isInt)).^2)).*sqrt(sum(abs(G.mean(isInt)).^2))); % coherence
            
            % Keep a few figures for all internals
            %FG(ci,:) = fg;
            %FGI(ci,:) = fg(isInt);
            %FGSEP(ci,:) = fgsep;
            %FGSEPI(ci,:) = fgsep(isInt);
            
            
            %% Combine internals for "average" phase diff and error
            switch cfg.intCombMethod
                case 'ampWeightSum'
                    % Method 1 - vector sum (i.e. amplitude weighted)
                    
                    PD(ci) = angle(I(ci));
                    Imd(ci) = mean(rangeCoarse(isInt).*abs(fg(isInt)))./mean(abs(fg(isInt))); % effective mean internal depth (amplitude weighted)
                    
                    % Now two methods for errors on internals
                    switch cfg.errorMode
                        case 'inBin'
                            % Calculate phase error by combining the errors from each bin (i.e. burst stats)
                            Ie(ci) = sqrt(sum(fgse(isInt).^2)); % total error of I from sum of noise within bins
                            PDe(ci) = Ie(ci)./(sqrt(2)*abs(I(ci))); % phase error of I
                        case 'btwnBin'
                            % phase-difference variance between internals
                            % done through the ESA equation for phase variance based on coherence.
                            % Calculate phase variance across the internals from coherence
                            gamma = Ic(ci);
                            In = sum(isInt); % number of internals
                            PDe(ci) = (1/sqrt(2*In))*sqrt(1-abs(gamma)^2)./abs(gamma);% http://www.esa.int/esapub/tm/tm19/TM-19_ptC.pdf eq. 1.18
                    end
                    
                case 'errorWeightFit'
                    % Method 2 - estimate by direct error weighted fitting
                    x = transpose(real(fg(isInt)));
                    y = transpose(imag(fg(isInt)));
                    er = transpose(fgse(isInt));
                    if all(er==0)
                        er = ones(size(er));
                    end
                    %[gfgm,gfgsem] = menke_fit(x,y,er); %
                    %[gfgl,gfgsel] = lscov(x,y,1./er.^2);
                    %[gfg,stats] = robustfit(x,y,@(r)1./(r*er),1,'CONST') % this uses iteration
                    [gfg,stats] = robustfit(x,y,@(r)1./er.^2,1,'off'); % weighted ols
                    gfgse = stats.se; % 95% ci
                    
                    % Note ideally we'd use total least squares, and this
                    % technique would falll over if PD = pi/2, but in
                    % reality PD is usually close to 0;
                    PD(ci) = atan(gfg); % mfg is slope of fg, i.e. d(imag(fg))./d(real(fg))
                    mpred = x + 1i*gfg*x;
                    mpredse = gfgse*x;
                    %resid = fg(isInt)-mpred;
                    %PDe(ci) = std(angle(fg(isInt))-PD(ci)); % emperical phase std (not nice having std of angle... but ok for small PD)
                    PDe(ci) = atan(gfgse); % phase std (as calculated by menke fit)
                    Imd(ci) = mean(rangeCoarse(isInt)); % this not perfect - as actual range is weighted in menke fit...
                    
                case 'meanphasediff'
                    % Method 3 - ignore amplitude
                    % this is for the case where we have deterministic system
                    % noise causing phase difference offsets which are not
                    % related to signal strength. Here we can just 
                    pd = angle(fg(isInt))';
                    %PD(ci) = mean(pd);
                    PD(ci) = median(pd);
                    PDe(ci) = std(pd)./sqrt(length(pd)-1); % standard error
            end
            
            % Calculate noise floor
            %Fnm(ci) = mean(F.ste(isInt)); % mean noise
            %disp(['F noise mean: ' num2str(dB(Fnm(ci)))])
            %Fpnp = Fnm(ci)./(sqrt(2)*abs(F.mean(isInt))); % phase noise prediction assuming white noise in RF
            %Gnm(ci) = mean(G.ste(isInt)); % mean noise
            %Gpnp = Gnm(ci)./(sqrt(2)*abs(G.mean(isInt))); % phase noise prediction assuming white noise in RF
            % note there appears to be a minimum phase noise of 1e-3 rad...
            % wny - this should eb much lower?
            
            % Stop to investigate noise
            if cfg.doPlotTest && any(ci==cfg.plotCI) % (ci==10 || ci==49) % cfg.doPlotTest  %ci==N-1     || cfg.doPlotTest %
            %if cfg.doPlotTest && (ci==346 || ci==347) % cfg.doPlotTest  %ci==N-1     || cfg.doPlotTest %
                title_txt = [int2str(ci) ') F: ' F.filename ' b:' int2str(F.burst) ' t:' datestr(F.time)];
                
                % Amplitude profile and phase diff at internals
                figure
                ax(1) = subplottight(3,1,1);
                plot(rangeCoarse,dB(abs(F.mean)),'r')
                hold on
                plot(rangeCoarse,dB(abs(G.mean)),'b')
                %plot(rangeCoarse,dB(abs(fg)),'g')
                legend('f','g')
                set(gca,'XTickLabel',[])
                title(title_txt,'interpreter','none')
                
                ax(2) = subplottight(3,1,2);
                %plot(rangeCoarse(isInt),angle(fg(isInt)),'k.')
                [h3,h4] = erbar(rangeCoarse(isInt),PD(ci),-PDe(ci),PDe(ci),'g','k');
                hold on
                [h1,h2] = erbar(rangeCoarse(isInt),angle(fg(isInt)),fgsep(isInt),-fgsep(isInt),'r','b');
                
                ylim([-0.3 0.3])
                ylabel('phase difference')
                box on
                
                % Plot phase
                ax(3) = subplottight(3,1,3);
                plot(rangeCoarse,angle(F.mean),'r')
                hold on
                plot(rangeCoarse,angle(G.mean),'b')
                xlabel('range (m)')
                ylabel('phase')
                
                linkaxes(ax,'x')
                clear ax
                
                figure, 
                erbar(real(fg(isInt)),imag(fg(isInt)),fgse(isInt));
                hold on
                erbar(real(mpred),imag(mpred),-mpredse,mpredse,'g','k');
                plot(mpred,'g.')
                
                if 0
                    % scatter
                    figure
                    ax(1) = subplottight(2,1,1);
                    plot(dB(abs(F.mean(isInt))),F.stePhase(isInt),'b.','markersize',15);
                    hold on
                    plot(dB(abs(G.mean(isInt))),G.stePhase(isInt),'r.','markersize',15);
                    colormap jet
                    grid on
                    %hold
                    ylabel('phase standard error')
                    %ch = colorbar('East');
                    %ylabel(ch,'phase difference (rad)')
                    set(gca,'XTickLabel',[])
                    %ch.Label.string = 'amplitude (dB)';
                    
                    % overlay prediction of noise from this noise floor
                    hold on
                    %plot(dB(fliplr(sort(abs(F.mean(isInt))))),sort(Fpnp),'b'); % predicted phase noise for F
                    %plot(dB(fliplr(sort(abs(G.mean(isInt))))),sort(Gpnp),'r'); % predicted phase noise for F
                    %plot(dB(abs(G.mean(isInt))),Gpnp,'r.'); % predicted phase noise for G
                    xlabel('amp')
                    title(title_txt)
                    
                    % Plot phase difference
                    ax(2) = subplottight(2,1,2);
                    erbar(dB(sqrt(abs(fg(isInt)))),angle(fg(isInt)),fgsep(isInt),-fgsep(isInt),'r','b');
                    xlabel('amp')
                    ylabel('phase difference')
                    linkaxes(ax,'x')
                    
                    keyboard
                end
            end
        end
    end
end

% Strain
Br = transpose(rangeCoarse(bedInd)); % bed range (through time
IBr = Br-Imd; % internal - bed range
drs = IBr.*dt*(cfg.vsr(1)/daysPerYear); % range change from strain
drse = IBr.*dt*(cfg.vsr(2)/daysPerYear); % range change error from strain

% Range change
dr = F.lambdac*PD./(4*pi);
dre = F.lambdac*PDe./(4*pi);
tr = -daysPerYear*(dr)./dt; % thinning rate
mr = -daysPerYear*(dr-drs)./dt; % melr rate (thinning corrected for strain)
mre = -daysPerYear*sqrt(dre.^2 + drse.^2)./dt; % Combining radar error with strain error

% r = cumsum(dr); % range

% % Smoothing...
% % Using singular spectrum analysis
% mrdt = detrend(mr);
% mrt = mr-mrdt; % trend
% n = round(length(mr)/2);
% mrssa = ssa(mrdt,n,'noplot','silent');
% n2 = round(0.8*n);
% mrsdt = ssarecon(mrssa,n2:n,'noplot','silent');
% mrs_ssa = mrsdt + mrt; % retrended to reconstruct smootghed meltrate time series

% Running mean
%mrs3 = runningmean(mr,3);
%mrs5 = runningmean(mr,5);

% cfg.lambdac = vdat.lambdac;

% %% Save results to file
[~,cfgfilename,~] = fileparts(cfgfile);
outfilename = ['fmcw_timeseries_' cfgfilename];
disp(['Saving data to: ' outfilename])
save(outfilename,'timeC','Ic','Imd','IBr','drs','drse','PD','PDe','rangeCoarse','mr','mre','tr','dr','dre','bedInd','cfg') % 'temp1','temp2'

% Calculate melt and plot
%fmcw_ts_clean(outfilename)

fmcw_ts_plot(outfilename)

function X = getSpec(filename,burst,cfg)

[vdat] = fmcw_load(filename,burst); % load data

switch cfg.attSetKeep
    case 'min'
        % Now keep only chirps with the first attenuator setting
        vdats = fmcw_burst_split_by_att(vdat);
        attsum = vdat.Attenuator_1 + vdat.Attenuator_2;
        [~,minatti] = min(attsum);
        %disp(['min att num is: ' int2str(minatti)])
        vdat = vdats(minatti); % Just use the min attenuation settting for simplicity (for now)
    case 'all'
        % Do northing!
end

% Keep only selected chirp
if ~isempty(cfg.chirpNo)
    vdat = fmcw_burst_subset(vdat,cfg.chirpNo);
end

%% Cleaning
% Frequency selection
if min(cfg.prepro.freqUseRange)>2e8 || max(cfg.prepro.freqUseRange)<4e8
    vdat = fmcw_cull_freq(vdat,cfg.prepro.freqUseRange);
end
    
% Drop first N chirps
if cfg.prepro.cullFirstChirps > 0
    chirpsToKeep = [cfg.prepro.cullFirstChirps+1:vdat.ChirpsInBurst];
    vdat = fmcw_burst_subset(vdat,chirpsToKeep);
end
    
% Kill outliers
%cfg.prepro.chirp.cullBad
%noisePowerLimit = 0.01;
%[vdat,~] = fmcw_cull_bad(vdat,noisePowerLimit,0);
    
% Keep clean chirps only
if ~isempty(cfg.prepro.keepCleanChirps)
    vdat = fmcw_keep_clean_chirps(vdat,cfg.prepro.keepCleanChirps,0);
end

% Process all chirps
%c = fmcw_range2(vdat,cfg.padfactor,cfg.bedSearchRange(2));
[c.rangeCoarse,~,~,c.specRaw] = fmcw_range(vdat,cfg.padfactor,cfg.bedSearchRange(2));

% Configure output
X.filename = filename;
X.rangeCoarse = c.rangeCoarse;
X.specRaw = c.specRaw;
X.time = vdat.TimeStamp;
X.temp1 = vdat.Temperature_1;
X.temp2 = vdat.Temperature_2;
X.lambdac = vdat.lambdac;
X.burst = burst;

function S = derampStats(A,bedInd)
% subtract phase from bed phase for each chirp and calculate burst stats

% Subtract the phase at depth from the bed phase (for each chirp)
specRB = conj(A.specRaw).*repmat(A.specRaw(:,bedInd),1,size(A.specRaw,2))./abs(mean(A.specRaw(:,bedInd))); 
% Take mean and get stats
S = fmcw_burst_stats(specRB); 
S.bedInd = bedInd;
% Now pass a few fields through that we need later in S
S.time = A.time;
S.filename = A.filename;
S.burst = A.burst;
S.lambdac = A.lambdac;
S.temp1 = A.temp1;
S.temp2 = A.temp2;