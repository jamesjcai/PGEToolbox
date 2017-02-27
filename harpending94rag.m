function [r]=harpending94rag(hap,showhist)
%HARPENDING94RAG - 
%
% [r]=harpending94rag(hap,showhist)
%
% Mismatch is describing pairwise differences between haplotypes but
% raggedness is the variation around the curve. An empirical mismatch
% distribution that does not deviate from a unimodal distribution of
% pairwise differences among haplotypes and has a smooth distribution
% (Harpending, 1994) suggests recent population expansion (Slatkin and
% Hudson, 1991; Rogers and Harpending, 1992). In other words, a mismatch
% distribution, P > 0.05 means you can’t reject the null hypothesis of
% population expansion.
%
% See also: MISMCH

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<2
    showhist=false;
end
[r]=raggedness(hap,showhist);

