function [ldinfo] = snp_ldpairh(haplodata)
%SNP_LDPAIRH - calculate pairwise LD from haplotype data
%
% [ldinfo] = snp_ldpairh(haplodata)
%
%
% SEE ALSO: EMLDRUN, SNP_LDPLOT, SNP_LDPAIR

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


[D_raw,D_prime,R2]=linkdisequ(haplodata,1);

ldinfo=struct;
ldinfo.d=D_raw;
ldinfo.dprime=D_prime;
ldinfo.r2=R2;