function [d] = snp_interpopdiv(geno1, geno2)
%SNP_INTERPOPDIV
%
%http://www.genetics.org/cgi/reprint/170/3/1181
%
% See also: snp_fst(geno1,geno2,'hughes')

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

[maf1, allemajor1, alleminor1] = snp_maf(geno1);
[maf2, allemajor2, alleminor2] = snp_maf(geno2);

idx = find(allemajor1 ~= allemajor2);
maf2(idx) = 1 - maf2(idx);

d = 1 - (sqrt(maf1.*maf2) + sqrt((1 - maf1).*(1 - maf2)));
