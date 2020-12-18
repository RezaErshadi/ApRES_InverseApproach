function fmcw_ts_plot(filename)

% Plot fmcw_timeseries as processed by fmcw_ts_range3

if nargin == 0
    filename = uigetfile('fmcw_timeseries_*.mat','Select a processed times series file');
    if filename==0
        return
    end
    cmd = [mfilename '(''' filename ''')'];
    disp('To repeat this processing paste the following command from the clipboard')
    disp(cmd)
    clipboard('copy',cmd)
end
load(filename)

%% Plot

% Internal and bed phase
figure
bx(1) = subplottight(3,1,1);
%h1 = plot(timeC,PD,'col','b');
%[hb,hp] = erbar(timeC,PD,-PDe,PDe,'k','k');
h = ershade(timeC,PD,-PDe,PDe);
hold on
h = plot(timeC,PD,'col','k');
grid on
ylabel('inter-shot phase difference (rad)')
legend(h,'bed-internal')
%legend('internal','bed','bed-internal')
xtimelab('tg')
% bx(3) = subplot(4,1,3);
% plot(timeC,dr);
% ylabel('range change (m)')
% xtimelab

% bx(2) = subplottight(3,1,2);
% erbar(timeC,mr,mre,-mre,'r','b'); % [hb,hp] =
% h = ershade(timeC,PD,-PDe,PDe);
% hold on
% h1 = plot(timeC,mr,'col',[0.6 0.6 0.6]);

% Plot smoothed melt rates
%hold on
%plot(timeC,mrs_ssa,'k')
%h2 = plot(timeC,mrs3,'g');
%h3 = plot(timeC,mrs5,'c');
%legend([h1 h2 h3],'raw','3-point mean','5-point mean')
%grid on
%ylabel('meltrate (m/year)')
%xtimelab('tg')

% melt rate
bx(2) = subplottight(3,1,2);
%h = ershade(timeC,mr,-3*mre,3*mre,[0.8 0.8 0.8]);
hold on
%h = ershade(timeC,mr,-2*mre,2*mre,[0.6 0.6 0.6]);
h = ershade(timeC,mr,-mre,mre,[0.3 0.3 0.3]);
h0 = plot(timeC,tr,'col','r');
h1 = plot(timeC,mr,'col','k');
legend([h0 h1],'thinning','melt')
%h2 = plot(timeC,mrs3,'g');
%h3 = plot(timeC,mrs5,'c');
%legend([h0 h1 h2 h3],'thinning','melt','melt 3-point mean','5-point mean')
grid on
ylabel('meltrate (m/year)')
xtimelab('tg')





%
% % Errors
% % Time series of phase errors and bed index
% figure
% ax(1) = subplot(2,1,1);
% plot(timeC,Ipe,'b')
% %legend('internals','bed')
% xtimelab('tg')
% ylabel('phase error')
%
% ax(2) = subplot(2,1,2);
% plot(timeC,bedInd,'r')
% ylabel('bed index')
% xtimelab
% linkaxes(ax,'x')

% % Time series of range and range error
% figure
% ax(1) = subplot(2,1,1);
% plot(timeC,r,'-.')
% ylabel('Cumulative range change (m)')
%
% ax(2) = subplot(2,1,2);
% plot(timeC,dre,'r')
% ylabel('range error')
% linkaxes(ax,'x')

% Timeseries of coherence
bx(3) = subplottight(3,1,3);
plot(timeC,abs(Ic),'.')
ylabel('coherence')
xtimelab

linkaxes(bx,'x')


% % Timeseries of phase-difference for each internal
% figure
% % add errorbars
% numInt = sum(isInt);
% for ii = 1:numInt
%     hold on
%     h(ii) = plot(timeC,angle(FGI(:,ii)));
%     %col = get(h(ii),'col');
%     legtxt(ii) = {[num2str(intDepth(ii)) ' m']};
%     %erbar(timeC,angle(FGI(:,ii)),-FGSEPI(:,ii),FGSEPI(:,ii),'col',col);
% end
% xtimelab
% ylabel('phase difference')
% legend(h,legtxt)

% % Noise floor timeseries
% figure
% plot(time(1:end-1),dB(Fnm),'b')
% hold on
% %plot(timeC,dB(Gnm),'r')
% xtimelab
% ylabel('noise floor')
