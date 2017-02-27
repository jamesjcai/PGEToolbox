function [p] = snp_genoprct(geno,dim)
%SNP_GENOPRCT - returns percentage of valid genotype data
%  Syntax: [p] = snp_genoprct(geno)
%
% i.e., Percent successfully genotyped samples
%
%  Syntax: [p] = snp_genoprct(geno,dim)
%                dim=1, marker-wise; dim=2, sample-wise


% Other:
%Percent successfully genotyped samples (Sample genotype success rate)
%Average genotyping success rate
%Duplicate sample error rate
%Non-Mendelian inheritance error rates (errors not consistent with normal
%transmission of chromosomes in family members)


% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if nargin<2,
    dim=1;
end

[n,m]=size(geno);
if(mod(m,2)>0), error(''); end

switch dim
    case 1
        locnum=m/2;
        p=zeros(1,locnum);
        parfor k=1:locnum
            locgeno=geno(:,[k*2-1 k*2]);
            locgeno=locgeno(:);
            p(k)=i_genopct(locgeno);
        end
    case 2
        idvnum=n;
        p=zeros(1,idvnum);
        parfor k=1:idvnum
            locgeno=geno(k,:);
            p(k)=i_genopct(locgeno);
        end
end


function [p] = i_genopct(locgeno)
    n=length(locgeno);
    x=sum(locgeno==5|locgeno==0);
    p=(n-x)./n;



% function [p] = i_genopct_old(locgeno)
% 	[numHap,sizHap,seqHap] = counthaplotype(locgeno);
% 	[n,m]=size(seqHap);
% 	existing=0;
% 	missing=0;
% 	for (k=1:n),
% 		allele1=seqHap(k,1);
% 		allele2=seqHap(k,2);
% 	      if (allele1~=5 && allele2~=5),
%  		     existing=existing+sizHap(k);
% 	      else
% 	             missing=missing+sizHap(k);
% 	      end
% 	end
% 	p=existing/(existing+missing);


