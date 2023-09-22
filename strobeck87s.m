function [StrobeckS] = strobeck87s(aln)
%STROBECK87S - Strobeck's S statistic
% Syntax: [StrobeckS] = strobeck87s(aln)

%The Strobeck's S test statistic (Strobeck 1987; see also Fu 1997) is also
%based on the haplotype (gene) frequency distribution conditional the value
%of theta (Ewens 1972, equations 19-21). The S statistic gives the probability
%of obtaining a sample with equal or less number of haplotypes than the observed.

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

%test statistics (Strobeck 1987; Fu 1996, 1997) use the sampling distribution of
%K directly .  number of haplotypes (K)

[Fs, StrobeckS] = fu97Fs(aln);
