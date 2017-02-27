function [markn,smpln]=snp_marklen(geno)
%SNP_MARKLEN - Number of markers
%
%See also: SNP_SAMPLEN, SNPMARKBPLEN

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

[smpln,m2]=size(geno);
markn=m2/2;
