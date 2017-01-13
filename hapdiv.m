function [hd]=hapdiv(sizHap,unbiased)
%HAPDIV - Haplotype diversity/heterozygosity (hd)
%
%Returns unbiased hd (Nei Eq 8.4, p.178) and var(hd) (Nei Eq 8.12, p.180)
%also as described in Nei & Tajima (1981).
%if 'biasflag' is set >0, the biased version (biased by a factor of
%(n-1)/n) of the statistic is returned.
%
%Depaulis F, Veuille M.
%Neutrality tests based on the distribution of haplotypes under an infinite-site model.
%Mol Biol Evol. 1998 Dec;15(12):1788-90.
%PMID: 9917213

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if(nargin<2), unbiased=true; end
n=sum(sizHap);

%frqHap=sizHap/sum(sizHap);
%hd=1-sum((sizHap./m).^2);
hd=1-sum(sizHap.^2)./n./n;

if (unbiased),
    hd=hd/(1-1./n);
end

%hd=1-sum((ni./n).*((ni-1)./(n-1)));

