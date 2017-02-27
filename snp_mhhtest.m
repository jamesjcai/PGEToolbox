function [rmhh,rhh]=snp_mhhtest(geno1,geno2)
%SNP_MHHTEST - MHH Test, most frequent haploytpe homozygosity test
% Syntax: [rmhh,rhh]=snp_mhhtest(geno1,geno2)
%
%Kimura R, Fujimoto A, Tokunaga K, Ohashi J.
%A practical genome scan for population-specific strong selective sweeps that have reached fixation.
%PLoS ONE. 2007 Mar 14;2:e286.
%PMID: 17356696

%If a region showing MHH ?0.9 is unusually extended in the test population,
%rHH exhibits a small value. Therefore, rHH represents the extent of the
%population-specific decrease in haplotype diversity and thus can be an
%index for detecting a recent rapid increase in the frequency of an allele.
%Here, we can control the heterogeneity of recombination rates if the
%reference population can be regarded as neutral. However, rMHH and rHH
%would not show small values in the loci where the same allele is selected
%and fixed in both populations used in the comparison. In case different
%mutations in the same region were selected in the two populations, rMHH
%would show a small value, but rHH would not.

%rMHH<0.05 and rHH<0.3 corresponded to approximately 90% detection power
%(81.5%–91.5% for rMHH and 86.8%–94.5% for rHH) (Fig 2G and H). Therefore,
%we used these values as thresholds in subsequent analyses. Type I error
%rates for rMHH<0.05 and rHH<0.3 ranged from 0.24%–1.43% and 0.76%–4.62%,
%respectively.

% Simulation analyses suggested that rMHH has a higher ability to
% capture highly differentiated regions than maxFST does when the density
% of the typed SNPs is low

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

popid1='';
popid2='';
if nargin<1
       [marker,popid1,popid2]=selectMarker2Popcodes('EDAR',5,1);
       if isempty(marker)
           return;
       end

disp('Downloading GENODATA from population 1....')
       [geno1,mark1] = snp_downloadhapmap(marker,popid1);
       [geno1]=snp_breaktrio(geno1,mark1);

disp('Downloading GENODATA from population 2....')
       [geno2,mark2] = snp_downloadhapmap(marker,popid2);
       [geno2]=snp_breaktrio(geno2,mark2);

disp('Extracting common markers....')
x=intersect(mark1.pos,mark2.pos);
[a]=find(ismember(mark1.pos,x));
[b]=find(ismember(mark2.pos,x));
[geno1,mark1]=snp_pickmarker(geno1,mark1,a);
[geno2,mark2]=snp_pickmarker(geno2,mark2,b);

end



figure;
n=snp_samplen(geno1)./2; m=snp_marklen(geno1);
x=1; y=2; if n<m, x=2; y=1; end

subplot(x,y,1);
%figure;
snp_vgview(geno1); title(sprintf('Test Population %s\n',popid1));

subplot(x,y,2);
%figure;
snp_vgview(geno2); title(sprintf('Reference Population %s\n',popid2));

%subplot(3,1,3); plot(rand(1,snp_marklen(geno1))); xlim([1
%snp_marklen(geno1)]);

mark1.rsid

i_dispheader('rMHH and rHH Tests*')
fprintf('MHH: most frequent haplotype homozygosity\n');
fprintf(' HH: haplotype homozygosity\n\n');
%fprintf('rMHH: ratio of MHHs between test and reference populations\n');
%fprintf('rHH: ratio of HHs between test and reference populations\n');
fprintf('Test Population %s:\n ',popid1);
[mhh1,hh1]=snp_mhh(geno1);
fprintf('Reference Population %s:\n ',popid2);
[mhh2,hh2]=snp_mhh(geno2);
fprintf('------------------------------------------\n');

if mhh1>0
    rmhh=mhh2/mhh1;
    if rmhh<0.05
    fprintf('rMHH = %f (significant rMHH<0.05)\n',rmhh);
    else
    fprintf('rMHH = %f (not significant rMHH>=0.05)\n',rmhh);
    end
else
    fprintf('rMHH = NaN (not significant)\n');
end

if hh1>0
    rhh=hh2/hh1;
    if rhh<0.3
    fprintf(' rHH = %f (significant rHH<0.3)\n',rhh);
    else
    fprintf(' rHH = %f (not significant rHH>=0.3)\n',rhh);
    end
else
    fprintf(' rHH = NaN (not significant)\n');
end
i_dispfooter

disp('*REF: Kimura R, Fujimoto A, Tokunaga K & ')
disp('      Ohashi J., PLoS ONE. 2007 Mar 14;2:e286.')


