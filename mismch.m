function [Y] = mismch(aln)
%MISMCH - mismatch distribution
%
%distribution of pairwise sequence differences
%
%Count the number of site differences between each pair of sequences in a sample, and use the
%resulting counts to build a histogram.
%

%theta=1; a=0:25;           % theta=2Nu, a is the mutation difference
%f=(1/(theta+1))*((theta/(theta+1)).^a);        % f is fraction of
%expectation
% REF: Slatkin and Hudson (1993)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (isstruct(aln)), seq = aln.seq;
else seq = aln;
end

%[n,m] = size(seq);
%for i=1:n-1
%for j=i+1:n
%	D(i,j)=sum(seq(i,:)~=seq(j,:));
%end
%end

[n, p] = size(seq);
k = 1;
Y = zeros(1, n*(n - 1)./2);
for i = 1:n - 1
    dsq = zeros(n-i, 1);
    for q = 1:p
        dsq = dsq + (seq(i, q) ~= seq((i + 1):n, q));
    end
    Y(k:(k + n - i - 1)) = dsq;
    k = k + (n - i);
end
%squareform(Y)


%http://www.genetics.org/cgi/reprint/129/2/555.pdf
%We conclude that plotting frequency distribution of
%pairwise differences in sequence or equivalently pairwise
%divergence times of genes sampled provides an
%indication of the structure of the phylogenetic tree
%representing the history of those genes. It is difficult,
%however, to compare such a graph with a geometric
%distribution and reject the null hypothesis that the
%genes sampled were from a randomly mating population
%of constant size.

%http://mbe.oxfordjournals.org/cgi/reprint/13/7/895.pdf
%Mismatch distributions are histograms showing the pattern of nucleotide (or restriction) site differences between
%pairs of individuals in a sample. They can be used to test hypotheses about the history of population size and
%subdivision (if selective neutrality is assumed) or about selection (if a constant population size is assumed).
