function snp_haphomozygosity(hap, centralmark)
%SNP_HAPHOMOZYGOSITY - A simple pairwise metric of haplotypic homozygosity
%
%"Consider two biallelic SNPs, denoted A and B, where we are particularly
% interested in allele X of SNP A. We define the following metric:
%    Hx = h_BX - h_B
% where h_BX is the homozygosity observed at SNP B when we consider only
% haplotypes carrying allele X of SNP A, and h_B is the homozygosity
% observed at SNP B when we consider all haplotypes. We are particularly
% interested in situations where the homozygosity of the partitioned
% haplotypes is unusually high or low compared with the result expected
% from the general population. Therefore, we calculated a calibrated
% partition homozygosity (Hx) by subtracting the general-population
% homozygosity at SNP B (h_B) from h_BX.
%
% A positive value of Hx implies that the homozygosity at locus B on
% haplotypes marked by allele X of SNP A is greater than would be expected,
% given the population allele frequency of SNP B. An important feature of
% this metric is that each allele of a single SNP is likely to receive a
% different Hx value, corresponding to the history of that allele."
%
%Reference:

%{
A window size of 50 kb centered on the target allele contains sufficient
markers (?20 SNPs) to demonstrate our signal clearly; however, alternative
window sizes could be employed.

Fry AE, Trafford CJ, Kimber MA, Chan MS, Rockett KA, Kwiatkowski DP. Haplotype
homozygosity and derived alleles in the human genome. Am J Hum Genet. 2006
Jun;78(6):1053-9. Epub 2006 Apr 5. PubMed PMID: 16685655; PubMed Central PMCID:
PMC1474085.


p=(f11+f21)/(f11+f12+f21+f22);
px=f11/(f11+f12);
hB=1-(2*p*(1-p));
hBX=1-(2*px*(1-px));
Hx=hBX-hB;
%}


% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $
