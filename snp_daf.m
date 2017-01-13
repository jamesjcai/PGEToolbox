function [p,devalle] = snp_daf(geno,ancalle,ishaplotype)
%SNP_DAF - Derived Allele Frequence of SNP
%[p,devalle] = snp_daf(geno,ancalle)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-12-27 00:25:00 -0600 (Fri, 27 Dec 2013) $
% $LastChangedRevision: 755 $
% $LastChangedBy: jcai $

if nargin<3,
    ishaplotype=0;
end
if nargin<2
    warning ('No ancestral allele provided. SNP_DAF computes MAF instead!')
    [p,~,alleminor] = snp_maf(geno,ishaplotype);
    devalle=alleminor;
    return;
end




switch ishaplotype

case 0
    m2=size(geno,2);
    if(mod(m2,2)>0)
        error('Wrong GENODATA!');
    end
    m=m2/2;
    if length(ancalle)~=m
        error('Wrong ANCALLE provided.')
    end
    p=zeros(1,m);
    devalle=zeros(1,m);

    for k=1:2:m2
          x=geno(:,[k k+1]);
          x=x(:);
          x(sum(x==5,2)>0,:)=[]; % remove
          y=sum(x==ancalle((k+1)/2))/length(x(:));
          p((k+1)/2)=1-y;
    end

   idx = ancalle<1 | ancalle>4;
   p(idx)=nan;

case 1
    
    % now treat geno as Haplotype data
    m=size(geno,2);
    if length(ancalle)~=m
        error('Wrong ANCALLE provided.')
    end

    p=zeros(1,m);
    devalle=zeros(1,m);

    for k=1:m
          x=geno(:,k);
          x(x==5)=[];
          y=sum(x==ancalle(k))/length(x);
          p(k)=1-y;
    end

end



