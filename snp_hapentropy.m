function [s] = snp_hapentropy(hap, winsize, stepsize, pos)
%SNP_HAPENTROPY - haplotype entropy
%
%Syntax: s=snp_hapentropy(hap)
%        figure; snp_hapentropy(hap,winsize)
%
% References:
% Mizuno H, Atwal G, Wang H, Levine AJ, Vazquez A. Fine-scale detection of
% population-specific linkage disequilibrium using haplotype entropy in the human
% genome. BMC Genet. 2010 Apr 23;11:27. PubMed PMID: 20416085; PubMed Central
% PMCID: PMC2873552.
%
% Nothnagel M, Fürst R, Rohde K. Entropy as a measure for linkage disequilibrium
% over multilocus haplotype blocks. Hum Hered. 2002;54(4):186-98. PubMed PMID:
% 12771551.
%
% Atwal GS, Bond GL, Metsuyanim S, Papa M, Friedman E, Distelman-Menachem T, Ben
% Asher E, Lancet D, Ross DA, Sninsky J, White TJ, Levine AJ, Yarden R. Haplotype
% structure and selection of the MDM2 oncogene in humans. Proc Natl Acad Sci U S A.
% 2007 Mar 13;104(11):4524-9. Epub 2007 Mar 5. PubMed PMID: 17360557; PubMed
% Central PMCID: PMC1838634

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if nargin < 3
    stepsize = 1;
end

if nargin < 2
    [numHap, sizHap, seqHap] = counthaplotype(hap);
    p = sizHap ./ sum(sizHap);
    s = -sum(p.*log2(p));
else
    if isempty(winsize)
        winsize = 21;
    end
    m = size(hap, 2);
    if m <= winsize
        error('Window size too big')
    end
    s = [];
    px = [];
    for k = 1:stepsize:m - winsize + 1
        hap1 = hap(:, k:k+winsize-1);
        s1 = snp_hapentropy(hap1);
        s = [s, s1];
        p1 = mean(pos(k:k+winsize-1));
        px = [px, p1];
    end
    if nargout < 1
        if pos
            plot(px, s);
        else
            plot(1:length(s), s);
        end
        xlabel('Position')
        ylabel('Haplotype Entropy')
    end
end
