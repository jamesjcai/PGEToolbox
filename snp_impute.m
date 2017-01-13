function [geno]=snp_impute(geno)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

[~,m2]=size(geno);
if(mod(m2,2)>0), error('Wrong GENODATA!'); end
m=m2/2;

for (k=1:2:m2),
      g=geno(:,[k k+1]);
      idx=find(g==5);
      
      if ~isempty(idx) 
         x=g(g~=5);
         a=min(x); b=max(x);
         if a==b
             g(idx)=a;
         else
             p=sum(g==a)./length(x);
             g(idx)=i_randimpute(p,a,b,length(idx));
         end
      end
     geno(:,[k k+1])=g;
end


function c=i_randimpute(p,a,b,n)
rand('seed',5489);
c=zeros(n,1);
for k=1:n
if rand>=p
    c(k)=a;
else
    c(k)=b;
end
end
    

