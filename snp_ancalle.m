function [ancalle,ancp]=snp_ancalle(geno)
%SNP_ANCALLE - returns major alleles as ancestral alleles
%
% [ancalle,ancp]=snp_ancalle(geno)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

[n,m2]=size(geno);
if(mod(m2,2)>0) 
    error('Wrong GENODATA!'); 
end
m=m2/2;

ancp=zeros(1,m);
%devalle=zeros(1,m);
ancalle=zeros(1,m);

for (k=1:2:m2),
      x=geno(:,[k k+1]);
      x=x(:);
      x(find(sum((x==5),2)>0),:)=[]; % remove
      [a,b]=UniqueValues(x);
      [b,idx]=sort(b);
      a=a(idx);
      ancalle((k+1)/2)=a(end);
      ancp((k+1)/2)=b(end)./length(x);      
end

