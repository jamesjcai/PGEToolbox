function [rmin, D] = fourgamete(seq, showit)
%Four-gamete test

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if nargin < 2, showit = false; end
[rmin, D] = hudsonkaplan85rm(seq, showit);