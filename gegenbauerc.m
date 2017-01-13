function [c]=gegenbauerc(n,m,x)
%GegenbauerC[n, m, x] gives the n-th Gegenbauer polynomial in x for parameter m.

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

cx = gegenbauer_poly(n,m,x);
c=cx(end);
