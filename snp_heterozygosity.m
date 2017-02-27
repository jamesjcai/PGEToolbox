function [h,hv]=snp_heterozygosity(geno)
%SNP_HETEROZYGOSITY - SNP heterozygosity 
%  Syntax: [h,hv]=snp_heterozygosity(geno) 
%
%Same as: SNP_DIVERSITY, but not give unbiased estimation
%
%Define the frequencies of two alternative nucleotides at SNP j to be pj
%and qj =1 - pj. The estimator for heterozygosity is 2piqi[n/(n - 1)],
%where n is the sample size. The heterozygosity for a window is the simple
%arithmetic average of heterozygosities of the SNPs in that window. 
%
%REF: http://www.genome.org/cgi/content/full/15/11/1496

%These heterozygosities were calculated for each of the three populations
%and were denoted at HSij, referring to the heterozygosity in subpopulation
%i for SNP j. For the pooled sample across the populations, we calculate
%the average (across populations) of allele frequencies (pj_bar and qj_bar), and define
%the total heterozygosity as H_Tj=2pj_bar*qj_bar. FST can be thought of as the component of
%variance in allele frequency that is between populations, and Wright's 
%approximate formula is FST =(HT - HS)/HT. 

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

n=snp_samplen(geno);
p=snp_maf(geno);
q=1-p;

x=n/(n-1);
hv=2.*x.*p.*q;

% Calculate the local expected heterozygosity, or gene diversity, formally
% use: 1-sum(p1^2+p2^2+p3^2...)
%
% With two alleles it is easier to use 2pq than to use the
% summation format above.

h=nanmean(hv);

if (nargout<1)
i_dispheader('SNP Heterozygosity')
	fprintf (['Mean SNP heterozygosity = %f\n'],h);
i_dispfooter
end

