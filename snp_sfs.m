function snp_sfs(genodata)
%SNP_SFS - Site-Frequency Spectrum of SNPs

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

x=snp_maf(genodata);
isfolded=1;
histsfs(x,isfolded);
