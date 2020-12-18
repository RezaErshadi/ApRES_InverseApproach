function vdat = fmcw_keep_clean_chirps(vdat,nKeep,doPlot)

% Keep only cleanest chirps from a FMCW radar burst
%
% args: vdat = fmcw data structure
% n = number of cleanest chirps to keep
%
% Craig Stewart
% 2014/5/20
% Adapted from fmcw_cull_noisey

% Check input args
if nargin <3
    doPlot = 0;
end

% Check if there are more than one attenuator settings
if length(unique(vdat.chirpAtt)) > 1
    disp('Warning: Multiple attenuator settings present in burst - subset into attenuator settings before culling noisy chirps')
    vdat.processing = [vdat.processing {[mfilename ': cancelled - multiple attenuator settings in burst']}];
    return
end

% Check we're not killing all the data
nchirps = size(vdat.vif,1);
if nKeep >= nchirps
    disp([mfilename ' cancelled  - already at desired number of chirps'])
    vdat.processing = [vdat.processing {[mfilename ': cancelled - already at desired number of chirps']}];
    return
end

% Iteratively remove noisest shot (i.e. boot strap)
isGood = true(nchirps,1); % initially flag all as good
p = ones(nchirps,1);
while sum(isGood)>nKeep
    % Find noisy shots
    mv = mean(vdat.vif(isGood,:),1); % mean of all good chirps
    goodChirpList = find(isGood);
    for cn = 1:length(goodChirpList)
        p(goodChirpList(cn)) = rms(vdat.vif(goodChirpList(cn),:)-mv); % rms difference from mean for each chirp
    end
   
    % Badflag noisiest chirp
    [~,bci] = max(p(goodChirpList)); % chirp to remove
    isGood(goodChirpList(bci)) = false;
end

% Remove noisiest chirp
chirpsToKeep = find(isGood);
vdat = fmcw_burst_subset(vdat,chirpsToKeep);
vdat.processing = [vdat.processing {[mfilename ': kept cleanest ' int2str(nKeep) ' chirps: ' mat2str(vdat.chirpNum)]}];

if doPlot
    figure
    h(1) = bar(p,'r');
    hold on
    h(2) = bar(isGood.*p,'g');
    legend('discarded','kept')
    xlabel('chirp')
    ylabel('volts rms')
    title('chirp noise level')
end