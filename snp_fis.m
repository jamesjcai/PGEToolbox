function [f] = snp_fis(geno)
%SNP_FIS -
%
%[f]=snp_fis(geno)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

n = snp_samplen(geno);
hi = snp_obshet(geno);
[hs1, hs2] = snp_diversity(geno);
f = (hs2 - hi) / hs2;
