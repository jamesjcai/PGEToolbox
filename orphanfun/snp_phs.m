function [p,coresets]=snp_phs(hapldata,n1,n2,showit)
%SNP_PHS - pairwise haplotype-sharing score
%
%  Syntax: [p,coresets]=snp_phs(hapldata,n1,n2,showit)
%
% n1   -   Position of first core marker
% n2   -   Position of last core marker
% P        - PHS value
% CORESETS - core haplotype set (used to be determined ancestral state)
%
% see also: SNP_EHH
%
%ref: Toomajian et al (2006) PLoS Biology

% This score includes a function that controls for population structure:
% since pairs of individuals from the same population are more closely
% related than those from different populations, they're more likely to
% share long haplotypes and could bias the results.
%
%Citation: Toomajian C, Hu TT, Aranzana MJ, Lister C, Tang C, et al. (2006)
%A Nonparametric Test Reveals Selection for Rapid Flowering in the
%Arabidopsis Genome . PLoS Biol 4(5): e137.
%doi:10.1371/journal.pbio.0040137







