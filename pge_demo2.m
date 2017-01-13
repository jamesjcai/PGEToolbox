%%PGEToolbox DEMO - SNP Analysis
% Welcome to PGEToolbox.  This is a demonstration of
% PGEToolbox's functions for SNP analysis.
%
% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $



cd(fileparts(which('pge_demo2')))
% Current working directory has been changed.

%%
% First, let's just load an example HapMap genotype data from human BRCA2
% locus

filename=fullfile('seq_examples','brca2.hmp');
[genodata,markinfo] = snp_readhapmap(filename);


%%
% Now let's view the marker information in this collection of SNPs
%
% Notice PGEToolbox display rsid, allele, chromosome, strand and position.

snp_viewmark(markinfo)

%%
% Now we plot relative position of SNPs on chromosome, the height of bars
% indicates MAF of SNPs.

markinfo.maf=snp_maf(genodata);
snp_plotmarkpos(markinfo)

%%
%SNP_VGVIEW brings a visual genotype (VG) view presenting complete raw
%datasets of individuals' genotype data. The display and interpretation of
%large genotype data sets can be simplified by this a graphical display.

cla;
snp_vgview(genodata)

%%
% SNP_PREDHET computes predicted percentage of heterozygosity individuals.

h=snp_predhet(genodata);
h'

%%
% Tajima's D test for SNPs.

snp_tajima89d(genodata);

%%
% SNPWEBLINK produces website links for given SNP markers.

snp_weblink(markinfo)

%%
% SNP_FREQPIE produces a pie chart of allele and genotype frequencies among
% populations for a given SNP

snp_freqpie

%%
% Extended Haplotype Homozygosity of phased haplotype data
% Here core is 28-30. Position of first core marker is 28,
% position of last core marker is 30.

figure;
load('example_data/phasedhaplotype','haplodata');
haplodata=haplodata;
snp_ehh(haplodata,128,128);


%%
% iHS (Integrated Haplotype Score) is a statistic that has
% been developed to detect evidence of recent positive
% selection at a locus
% Here we load a pre-computed example result file and plot
% scatter plot of genomic position and iHS.

figure;
load('example_data/ihs_result_example', 'ihs', 'markinfo');
ihs=ihs;
plot(markinfo.pos,ihs,'.');
ylabel('|iHS|')
xlabel('Genomic Postion')
title('iHS Scatter Plot')
hline(2.0)

%%
% Thank you for viewing this introduction to PGEToolbox SNP funcations.
