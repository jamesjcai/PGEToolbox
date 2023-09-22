function [d] = hapdiv_test(aln, unbiased)
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


if (nargin < 2), unbiased = 1; end
if (isstruct(aln)), seq = aln.seq;
else seq = aln;
end
n = size(seq, 1);
[~, sizHap] = counthaplotype(aln);
d = hapdiv(sizHap, unbiased);

if (nargout < 1)
    i_dispheader('Haplotype Diversity Test (Depaulis and Veuille, 1998)')
    fprintf('Number of Haplotypes, hd: %d\n', numHap);
    fprintf('Haplotype (gene) diversity, Hd: %f\n', d);

    thew = thetaw(aln);
    HD = hapdiv_simu(n, 1000, thew, 0);
    %p=sum(d<HD)./10000;
    p = 2 .* sum(HD > abs(d)) ./ 10000;
    fprintf('Statistical significance:\n P = %f%s\n', p, sigtag(p));
    i_dispfooter
end

%Watterson’s test is conditioned on the number of haplotypes, whereas Depaulis and Veuille's H
%is conditioned on S (and n).
