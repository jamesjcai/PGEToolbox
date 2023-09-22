function [haplodata] = snp_geno2hap(genodata)
%SNP_GENO2HAP - converts GENODATA to pseudo-HAPLODATA
% [haplodata] = snp_geno2hap(genodata)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


%n=snp_marklen(genodata)*2;   % how many alleles of SNPs

%n=size(genodata,2);
%haplodata=[genodata(:,[1:2:n]); genodata(:,[2:2:n])];


[n, m2] = size(genodata);
n2 = n * 2;
m = m2 / 2;
haplodata = zeros(n2, m);

haplodata([1:2:n2], :) = genodata(:, [1:2:m2]);
haplodata([2:2:n2], :) = genodata(:, [2:2:m2]);
