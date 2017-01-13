function [genodata,markinfo] = snp_hap2geno(haplodata)
%SNP_HAP2GENO - converts HAPLODATA to GENODATA
% [genodata,markinfo] = snp_hap2geno(haplodata)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


%n=snp_marklen(genodata)*2;   % how many alleles of SNPs


[n2,m]=size(haplodata);

if mod(n2,2)>0
    n2=n2-1;
    warning('One gamete is excluded.');
end

m2=m*2;
n=n2/2;

genodata=zeros(n,m2);
genodata(:,1:2:m2)=haplodata(1:2:n2,:);
genodata(:,2:2:m2)=haplodata(2:2:n2,:);

if nargout>1
    markinfo.rsid=num2cellstr(1:m);
    markinfo.pos=1:m;
end
    
