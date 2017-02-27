function [geno] = snp_12geno(geno)
%SNP_12GENO - Simplify GENO coding convention into 1 and 2. 

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

n=snp_marklen(geno);
for k=1:n
    gen=geno(:,2*k-1:2*k);
    genx=gen(:);
    genx(genx>4)=[];
    geny=gen;

   
    if max(genx)==min(genx)
        geny(:)=1;
    else
    [a]=unique(genx);
    if length(a)==2
        if(sum(genx==a(1))>=sum(genx==a(2)))  
            geny(gen==a(1))=1;
            geny(gen==a(2))=2;
        else
            geny(gen==a(1))=2;
            geny(gen==a(2))=1;
        end
    else
        error('Not all SNPs are biallelic!')
    end
    end
    geny(gen>4)=0;
    %geno2=[geno2,geny];
    geno(:,2*k-1:2*k)=geny;
end

