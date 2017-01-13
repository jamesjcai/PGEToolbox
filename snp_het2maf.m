function [p]=snp_het2maf(h)
%SNP_HET2MAF

% H    -  average heterozygosity of an SNP
% p    -  minor allele frequency.
%
%REF: (Zhang and Li 2005)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if h>0.5
    error('average heterozygosity <=0.5');
end
p=(1-sqrt(1-2*h))./2;


%h=(1-(1-2*p)^2)/2