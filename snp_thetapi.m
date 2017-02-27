function [thepi]=snp_thetapi(geno,mark,persite)
%SNP_THETAPI - theta-pi from SNPs
% Syntax: [thepi]=snp_thetapi(geno,mark,persite)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<3, persite=0; end
if nargin<2, mark=[]; end

[smpln]=snp_samplen(geno);

[p] = snp_maf(geno);
%p=p(p>0);   % maybe should get rid of MAF=0
q=1-p;

x=nansum(2.*p.*q);
y=smpln/(smpln-1);

thepi=x*y;
%thepi2 = sum(smpln^2 .* p .* q) / (smpln*(smpln-1)/2);   % just the same
% as above

if (persite),
   [L]=snp_markbplen(mark);
   thepi=thepi/L;
end
