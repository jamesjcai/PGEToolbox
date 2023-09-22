function [h] = snp_haphom(hap)

% Ref:
% Nei M 1975. Molecular population genetics and evolution. North-Holland/American Elsevier.
% Sabatti C & Risch N 2002. Homozygosity and linkage disequilibrium. Genetics 160: 1707-
% Sabeti PC et al. 2002. Detecting recent positive selection in the human genome from haplotype structure. Nature 419: 832-837

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


[numHap, sizHap, seqHap] = counthaplotype(hap);
p = sizHap ./ sum(sizHap);
h = (sum(p.^2) - 1 / numHap) ./ (1 - 1 / numHap);
