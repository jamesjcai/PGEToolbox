function [hh]=haphom(sizHap)
%HAPHOM - Haplotype homozygosity (HH)
%An effective measure of linkage disequilibrium (LD) for more than
%2 markers.
%
%See also: SNP_EHH, HAPDIV

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

%hh = [sum(pi.^2) - 1/n] / (1 - 1/n]

n=sum(sizHap(:));
p=sizHap./n;
hh=(sum(p.^2)-1/n)/(1-1/n);
