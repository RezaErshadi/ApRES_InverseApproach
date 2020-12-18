function y = runningmean(x,n)

% Filters the columns of x with running mean length n

% INPUTS:
% x: the signal for which we desire the running mean (one per column)
% n: size of the sliding window in which to compute the mean
%
% Craig Stewart
% 2013

if size(x,1)==1
    x = x';
    dotrans = 1;
end

coeffs = ones(1,n)./n;
for ii = 1:size(x,2)
    y(:,ii) = conv(x(:,ii),coeffs,'same');
end

if exist('dotrans','var')
    y = y';
end