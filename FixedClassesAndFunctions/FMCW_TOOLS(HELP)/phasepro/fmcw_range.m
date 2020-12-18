function [Rcoarse,Rfine,spec_cor,spec] = fmcw_range(vdat,p,maxrange,winfun)

% [Rcoarse,Rfine,spec_cor,spec] = fmcw_range(vdat,p,maxrange,winfun)
%
% Phase sensitive processing of FMCW radar data based on Brennan et al. 2013
%
% Based on Paul's scripts but following the nomenclature of:
% "Phase-sensitive FMCW radar imaging system for high precision Antarctic
% ice shelf profile monitoring"
% Brennan, Lok, Nicholls and Corr, 2013
%
% Summary: converts raw FMCW radar voltages into a range for 
%
% input args: 
% vdat = structure containing shot metadata
% vdat.vif = data matrix of vltages as digitised with one chirp per row size(nchirps,nsamples)
% p = pad factor (i.e. level of interpolation to use during fft)
% maxrange = maximum range to crop output to
% winfun = window function handle (defaults to blackman)
%
% outputs:
% Rcoarse = range to bin centres (m)
% Rfine = range to reflector from bin centre (m)
% spec_cor = spectrum corrected. positive frequency half of spectrum with
% ref phase subtracted. This is the complex signal which can be used for
% cross-correlating two shot segements.

% Craig Stewart
% 2013 April 24
% Modified frequencies 10 April 2014

if nargin < 3
    maxrange = 2000; %m (range to crop output to)
end
if nargin < 4
    winfun = @blackman; % default to blackman window
end

% Ensure odd number of samples per chirp (to get phase centering right)
vdat = fmcw_burst_make_length_odd(vdat);

% Extract variables from structure to make it readable
fs = vdat.fs;
T = vdat.T;
B = vdat.B;
K = vdat.K;
ci = vdat.ci;
fc = vdat.fc;
lambdac = vdat.lambdac;

% Processing settings
N = size(vdat.vif,2);
if mod(N,2)
    xn = 0.5*(N-1); % timeshift amount (steps) prior to fft to get fft measuring phase at t=T/2
    %shiftFrac = 0.45;
    %xn = round(shiftFrac*(N-1)); % timeshift amount (steps) prior to fft to get fft measuring phase at t=T/2
    %disp(['warning - non standard shift offset ' num2str(shiftFrac)])
else
    disp('Warning: even number of sample in record, not possible to correctly offset to phase centre')
    xn = 0.5*(N);
end
[nchirps,N] = size(vdat.vif);

% Measure the sampled IF signal: FFT to measure frequency and phase of IF
%deltaf = 1/(T*p); % frequency step of FFT
%f = [0:deltaf:fs/2-deltaf]; % frequencies measured by the fft - changed 16 April 2014, was %f = [0:deltaf:fs/2]; 
%Rcoarse = f*ci*T/(2*B); % Range at the centre of each range bin: eq 14 (rearranged) (p is accounted for inf)
%Rcoarse = [0:1/p:T*fs/2-1/p]*ci/(2*B); % Range at the centre of each range bin: eq 14 (rearranged) (p is accounted for inf)

nf = round((p*N)/2 - 0.5); % number of frequencies to recover
%nf = length(f); % changed from above 2014/5/22
%nf = length(Rcoarse); 
win = window(winfun,N); %chebwin(N);  %rectwin(N); %

%% Loop through for each shot in burst
% Calculate phase of each range bin centre for correction
n = (0:nf - 1)';
Rcoarse = transpose(n*ci/(2*B*p));
phiref = 2*pi*fc*n./(B.*p) - (K*n.^2)/(2*B.^2*p.^2); % eq 17: phase for each range bin centre (measured at t=T/2), given that tau = n/(B*p)
[spec,spec_cor] = deal(zeros(nchirps,nf)); % preallocate
for ii = 1:nchirps
    vif = vdat.vif(ii,:);
    vif = vif-mean(vif); % de-mean
    vif = win.*vif.'; % windowed
    %vif = [vif; zeros((p-1)*N,1)]; % zero padded to length p*N
    vifpad = zeros(p*N,1);
    vifpad(1:length(vif)) = vif;
    vifpad = circshift(vifpad,-xn); % signal time shifted so phase centre at start
    %plot(vifpad), keyboard
    fftvif = (sqrt(2*p)/length(vifpad)).*fft(vifpad); % fft and scale for padding 
    fftvif = fftvif./rms(win); % scale for window
    
    % Compensate phase for ADC delay of deltaT (by adding linear phase ramp)
    doSync = 0;
    if doSync
        if ii == 1;
            disp('Using sync tone for inter-chirp phase synchronisation')
        end
        if vdat.FileFormat == 2; % data from RMB-A prototype
            %Search for sync tone +-500 Hz of nominal
            fsync = 1e4;	% Frequency of sync tone
            [~,b1] = max(abs(fftvif(round((fsync-500)*p*T):round((fsync+500)*p*T)))); 
            b1 = b1 + round((fsync-500)*p*T) - 1; % index of sync tone in fftvif
            
            if 0 %ii == 1
                figure
                subplot(2,1,1), plot(dB(abs(fftvif)))
                hold on, plot(b1,dB(abs(fftvif(b1))),'ro')
                xlabel('sample num'), ylabel('magnitude')
                
                subplot(2,1,2), plot(angle(fftvif))
                hold on, plot(b1,angle(fftvif(b1)),'ro')
                xlabel('sample num'), ylabel('phase')
            end
            
            deltaT = -angle(fftvif(b1))/(2*pi*fsync); % phase delay (seconds)
            
        elseif vdat.FileFormat == 3 || vdat.FileFormat == 4; % data from RMB-A and B?
            syncStartOffset = 9154; % after chirp end ind
            syncEndOffset = 20000; % after chirp end ind
            syncTone = vdat.v(vdat.Endind(ii)+syncStartOffset:vdat.Endind(ii)+syncEndOffset);
            %figure, plot(syncTone)
            syncTone = syncTone-mean(syncTone);
            syncTone = fftshift(blackman(length(syncTone)).*syncTone);
            fftsync = fft(syncTone);
            [~,mai] = max(abs(fftsync(1:length(fftsync)/2-0.5))); %index of peak signal
            fsync = 1e4;	% Frequency of sync tone
            deltaT = -angle(fftsync(mai))/(2*pi*fsync); % phase delay (seconds)
            
%             % test plot
%             figure
%             ax(1) = subplot(2,1,1); % amp
%             plot(abs(fftsync),'b')
%             hold on
%             plot(mai,abs(fftsync(mai)),'ro')
%             ylabel('amp')
%             
%             ax(2) = subplot(2,1,2); % phase
%             plot(angle(fftsync),'b')
%             hold on
%             plot(mai,angle(fftsync(mai)),'ro')
%             ylabel('phase'), linkaxes(ax,'x')
%             keyboard
            
        end
        
        disp(['deltaT = ' num2str(deltaT)])
        
        % Calculate and apply the offset as a phase ramp to the spectrum
        % using nomenclature from http://en.wikipedia.org/wiki/Discrete_Fourier_transform#Shift_theorem
        k = (0:length(fftvif)-1)'; % zero base index of frequency
        m = deltaT*fs; % time shift in (fractional) sample bins
        % (note: padding in fft does not affect original sampling freq!)
        Ns = length(fftvif); % samples in fft
        deltaphi = 2*pi*m*k/Ns; % phase correction vector for fftvif (rad)
        fftvif = fftvif.*exp(1i*deltaphi); % correct raw spec for phase delay of this chirp
    end
    
    spec(ii,:) = fftvif(1:nf); % positive frequency half of spectrum up to (nyquist minus deltaf)
    comp = exp(-1i*(phiref)); % unit phasor with conjugate of phiref phase
    spec_cor(ii,:) = comp.*fftvif(1:nf); % positive frequency half of spectrum with ref phase subtracted
end
%Rfine = lambdac*angle(spec_cor)/(4*pi); % Distance from centre of range bin to effective reflector: eq 15
%Rfine = angle(spec_cor)./((4*pi/lambdac) - (4*Rcoarse*K/ci^2)); % this is the full equation including the term generated by the last term in (13)
Rfine = fmcw_phase2range(angle(spec_cor),lambdac,repmat(Rcoarse,size(spec_cor,1),1),K,ci);
%R = Rcoarse + Rfine;

% Crop output variables to useful depth range only
n = find(Rcoarse<=maxrange,1,'last');
Rcoarse = Rcoarse(1:n);
%Rcoarse = repmat(Rcoarse,nchirps,1); % make output same size as Rfine for consistence
Rfine = Rfine(:,1:n);
spec = spec(:,1:n);
spec_cor = spec_cor(:,1:n);
