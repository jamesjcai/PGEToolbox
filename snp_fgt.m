function [rmin,D]=snp_fgt(haplodata,showit)
%SNP_FGT - Four-gamete test for SNP haplotypes
% REF: (Hudson and Kaplan, 1985)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<2, showit=0; end

if nargout<1
    hudsonkaplan85rm(haplodata,showit);
else
    [rmin,D]=hudsonkaplan85rm(haplodata,showit);
end

