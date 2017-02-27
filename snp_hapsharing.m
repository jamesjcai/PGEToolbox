function [data]=snp_hapsharing(hap,pos,refalle)
%SNP_PHS - pairwise haplotype-sharing score
%
%  Syntax: [data]=snp_hapsharing(hap,pos,refalle)
%
% pos       -   Position of first core marker
% refalle   -   Reference allele
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

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

s=hap(:,pos);
idx=s==refalle;

n1=sum(idx);
n0=length(idx);

hapref=hap(idx,:);

z1=0;
for i=1:n1-1
for j=i+1:n1
    z1=z1+hapsharing(hapref(i,:),hapref(j,:),pos);
end
end
z1=z1./((n1*(n1-1))/2);

z0=0;
for i=1:n0-1
for j=i+1:n0
    z0=z0+hapsharing(hap(i,:),hap(j,:),pos);
end
end
z0=z0./((n0*(n0-1))/2);

data=z1-z0;













