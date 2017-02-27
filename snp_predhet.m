function [p] = snp_predhet(geno)
%SNP_PREDHET - predicted percentage of heterozygosity individuals

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

[n,m2]=size(geno);
if(mod(m2,2)>0), error('Wrong GENODATA!'); end
m=m2/2;
p=zeros(1,m);

for (k=1:2:m2),
      z=geno(:,[k k+1]);
      z(find(sum((z==5),2)>0),:)=[]; % remove
      [a,b,c]=unique(z(:));
      x=sum(c==1);
      y=sum(c==2);
      %p(1,(k+1)/2) = 1-(x^2+y^2)/(x+y)^2;
      p(1,(k+1)/2) = (2*x*y)/(x+y)^2;
end
