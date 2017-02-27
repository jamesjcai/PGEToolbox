function [bp] = chrlen(id)
%CHRLEN - returns human chromosome length
%USAGE: [bp] = chrlen(id)
%

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-11-24 17:13:39 -0600 (Sun, 24 Nov 2013) $
% $LastChangedRevision: 754 $
% $LastChangedBy: jcai $

bp=chrlen_hg19(id);
