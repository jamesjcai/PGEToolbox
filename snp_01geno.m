function [geno] = snp_01geno(geno,ancalle)
%SNP_01GENO - Simplify GENO coding convention into 0 and 1. 

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-04-29 14:33:14 -0500 (Mon, 29 Apr 2013) $
% $LastChangedRevision: 534 $
% $LastChangedBy: jcai $


if any(geno(:)>4),
    error('Not deals with undetermined alleles.')
end
if nargin<2, ancalle=[]; end

if isempty(ancalle)
    noanc=1;
else
    noanc=2;
end


n=size(geno,2)/2;

c=0;
for k=1:n
    gen=geno(:,2*k-1:2*k);
    geny=gen;
    
    switch noanc
    
    case 1
        genx=gen(:);
        %genx(genx>4)=[];
        if max(genx)==min(genx)
            geny(:)=0;
        else
            [a]=unique(genx);
            if length(a)==2
                if(sum(genx==a(1))>=sum(genx==a(2)))  
                    geny(gen==a(1))=0;
                    geny(gen==a(2))=1;
                else
                    geny(gen==a(1))=1;
                    geny(gen==a(2))=0;
                end
            else
                warning('Not all SNPs are biallelic!');
                c=c+1;
            end
        end
    
    case 2
        geny(gen==ancalle(k))=0;
        geny(gen~=ancalle(k))=1;
        
    end       
    
    %geny(geny>4)=0;
    geno(:,2*k-1:2*k)=geny;
end

c

    
    

