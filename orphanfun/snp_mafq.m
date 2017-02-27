function [p] = snp_mafq(geno)
%SNP_MAFQ - Quick Function for MAF of SNP
% [p] = snp_mafq(geno)
%
% This funciton take <4 mins to get MAF for HapMap SNPs
% tic; maf=snp_mafq(YRIgenodata); toc;
% Elapsed time is 202.808894 seconds.

% Population Genetics & Evolution Toolbox, (C) 2007
% Author: James J. Cai
% Email: jamescai@stanford.edu
% Website: http://bioinformatics.org/pgetoolbox/
% Last revision: 2/23/2007

if isempty(geno), p=[]; return; end
[n,m2]=size(geno);

if(mod(m2,2)>0), error('Wrong GENODATA!'); end
m=m2/2;
p=zeros(1,m);

for (k=1:2:m2),
      x=geno(:,[k k+1]);
      x=x(:);
      x(x>4|x<1)=[]; % remove
      
      cidx=(k+1)/2;
      
      if isempty(x)
        p(cidx)=nan;
      else
        
        a=max(x); b=min(x);
        if a==b
            y=0;
        else
          y=sum(x==a)/length(x);
   %       sum(x==a)
   %       length(x)
          y=min([1-y,y]);
        end     
        p(cidx)=y;
      end      
end




