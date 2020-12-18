% Copyright (C) 2001 Laurent Mazet
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% usage: b = de2bi(d, [n, [p]])
%
% Convert a non-negative integer to bit vector.
%
% d: positive integer
% n: number of rows of the ouput bit vector (default is max. size).
% p: base of decomposition (default is 2).
%
% b : bit vector.

% 2001-02-02
%   initial release

function b = de2bi(d, n, p)

  if (nargin == 1)
    p = 2;
    n = floor ( log (max (max (d), 1)) ./ log (p) ) + 1;
  elseif (nargin == 2)
    p = 2;
  elseif (nargin ~= 3)
    error ('usage: b = de2bi (d [, n [, p]])');
  end

  d = d(:);
  if ( any (d < 0) || any (d ~= floor (d)) )
    error ('de2bi: only handles non-negative integers');
  end

  power = ones (length (d), 1) * (p .^ [0 : n-1] );
  d = d * ones (1, n);
  b = floor (rem (d, p*power) ./ power);

end
