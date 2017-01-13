function [L]=snp_markbplen(markinfo)
%SNP_MARKBPLEN - genomic length in base pair where SNPs identified
%
%See also: SNP_SAMPLEN, SNPMARKLEN

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

L=max(markinfo.pos)-min(markinfo.pos)+1;