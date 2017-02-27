function [theh]=snp_thetah(geno,mark,persite,p,showwarning)
%SNP_THETAH - Fay's theta_H from SNPs
%
% Syntax: [theh]=snp_thetah(geno,mark,persite,p,showwarning)
%
% p -  derived allele frequency, see SNP_DAF

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<5, showwarning=1; end
if nargin<4, p=[]; end
if nargin<3, persite=0; end
if nargin<2, mark=[]; end

[smpln]=snp_samplen(geno);


if isempty(p)
    [p] = snp_maf(geno);
    if showwarning
        disp ('WARNING: SNP_THETAH uses MAF as derived frequency.')
    end
end


theh=nansum(2.*p.*p)*(smpln/(smpln-1));


if (persite),
   [L]=snp_markbplen(mark);
   theh=theh/L;
end
