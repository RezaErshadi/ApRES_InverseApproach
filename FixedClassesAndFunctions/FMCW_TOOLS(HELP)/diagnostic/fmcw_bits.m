function fmcw_bits(vdat)

% fmcw_bits(vdat)
%
% Display digitised bits for troubleshooting ADC

v = vdat.vif(1,:)'; % first chirp
%v = vdat.vif'; % all chirps
v = v(:); % column
c = 2^16*v./2.5; % rescale to 16 bit int
b = de2bi(c,16); % binary vector

figure
ax(1) = subplottight(2,1,1);
imagesc(b')
colormap gray
ylabel('bit (le)')
title(['Burst: ' int2str(vdat.Burst)]);

ax(2) = subplottight(2,1,2);
plot(c)
ylabel('adc count')
xlabel('sample')

linkaxes(ax,'x')

s = sum(abs(diff(b))); % number of sign changes per channel
figure
bar(1:16,s)
title(['Burst: ' int2str(vdat.Burst)])
ylabel('num changes')
xlabel('bit (le)')


% figure
% for ii = 1:16;
%     ax(ii) = subplottight(4,4,ii);
%     plot(b(:,ii))
%     legend(['bit ' int2str(ii)])
% end
%linkaxes(ax,'x')
