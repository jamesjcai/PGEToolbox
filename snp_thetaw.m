function [thew]=snp_thetaw(geno,mark,persite)
%SNP_THETAW - theta-W from SNPs
% Syntax: [thew]=snp_thetaw(geno,mark,persite)

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

[Sn]=snp_segsites(geno);

nx=1:(smpln-1);
a1 = sum(1./nx);
thew=Sn/a1;

if (persite),
   [L]=snp_markbplen(mark);
   thew=thew/L;
end


%if (nargout>1),
%    a2 = sum(1./nx.^2);
%    thewvar = thew./a1+(a2*thew^2)./a1^2;
%if persite
%    thewvar = m*thew./a1+(a2*(thew*m)^2)./(m*a1)^2;
%end
%end
