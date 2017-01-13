function [hapsim]=snp_hapsimilarity(haplodata,winsize,showit)
%SNP_HAPSIMILARITY - Haplosimilarity
% Syntax: [hapsim]=snp_hapsimilarity(haplodata,winsize)
%
% input: 
%         winsize = sliding window size (number of SNPs)
%         
%Reference:
% Hanchard NA, Rockett KA, Spencer C, Coop G, Pinder M, Jallow M, Kimber M, McVean G, Mott R, Kwiatkowski DP.	
% Screening for recently selected alleles by analysis of human haplotype similarity.
% Am J Hum Genet. 2006 Jan;78(1):153-9.

%PMID: 16385459
%see also: http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=1474085

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (nargin<3), showit=false; end
if nargin<2||isempty(winsize), winsize=10; end

[m]=size(haplodata,2);
if m<winsize
    error('Too few SNPs or too big window size.')
end

hapsim=ones(1,m-winsize+1);
for k=1:m-winsize+1
 hap=haplodata(:,k:k+winsize-1);
 [numHap,sizHap,seqHap]=counthaplotype(hap(:,1));   % first SNP
 
 if numHap>1
     xallele=seqHap(2,1);
     [numHap,sizHap,seqHap]=counthaplotype(hap);
     idx=find(seqHap(:,1)==xallele);
     p=sizHap(idx)./sum(sizHap(idx));
     hapsim(k)=sum(p.^2);
 end
end

if (nargout<1 || showit)
    i_dispheader('Haplosimilarity (Hanchard et al. 2006)')
    fprintf('Window size (# of SNPs): %d\n', winsize);
    fprintf('Haplosimilarity score (unnormalized): %f\n', sum(hapsim));
    figure;
    subplot(2,1,1)
    plot(hapsim)
    xlim([1,length(hapsim)])
    ylabel('Haplosimilarity score (unnormalized)')
    xlabel('Markers')
    i_dispfooter
end

