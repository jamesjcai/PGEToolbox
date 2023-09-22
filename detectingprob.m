function [u] = detectingprob(p, n)
%DETECTINGPROB - Probability of detecting both alleles in SNPs discovery.
%
%[u]=detectingprob(p,n)
%
% the probability that we will observe a mutation in a sample of n
% sequences  [http://www.genetics.org/cgi/content/full/162/4/2017]
%
%
%For any given SNP with known allele frequency P,
%the probability of detecting both alleles in an original sample
%of N chromosomes is u=1-p.^n-(1-p).^n.
%
%Example:
%
%p=0:0.01:0.5;
%n=[2 8 16 24 48 96];
%figure;
%hold on
%for k=1:numel(n)
%plot(p,detectingprob(p,n(k)));
%end
%xlabel('Minor Allele Frequency')
%ylabel('Fraction of SNPs Discovered')
%hold off
%
%See Also: FIXPROB, polymorphicpdf

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

u = 1 - p.^n - (1 - p).^n;
