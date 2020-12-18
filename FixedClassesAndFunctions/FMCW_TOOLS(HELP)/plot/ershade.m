function h = ershade(x,y,eh,el,patchcol,varargin)

% [hb,hp] = ershade(x,y,eh,el,patchcol,varargin)
%
% simple error shading, inputs must be vectors
%
% Craig Stewart
% 2014/11/16

if nargin<4
    el = -eh; % symetrical errors
end
if nargin<5
    patchcol = [0.9 0.9 0.9];
end

x = reshape(x,1,numel(x));
y = reshape(y,1,numel(y));
eh = reshape(eh,1,numel(eh));
el = reshape(el,1,numel(el));
if length(eh) == 1
    eh = repmat(eh,size(x));
end
if length(el) == 1
    el = repmat(el,size(x));
end

% Plot
h = patch([x fliplr(x)],[y+eh fliplr(y+el)],patchcol,'edgeColor','none',varargin{:}); % 