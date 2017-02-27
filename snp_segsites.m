function [Sn]=snp_segsites(geno)
%SNP_SEGSITES - number of polymorphis sites of SNPs

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

[p] = snp_maf(geno);
% Only count sites where there is actual diversity in the population
% Sn = sum((p~=0).*(p~=1));
Sn = sum(p>0 & p<1);