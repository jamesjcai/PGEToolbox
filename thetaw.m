function [t, tvar] = thetaw(aln, fromS, persite)
%THETAW - Watterson's theta from the number of segregating sites
%
% Syntax: [t,tvar] = thetaw(aln,fromS,persite)
%
%Watterson, G. (1975) On the number of segragating sites in genetical models without recombination. Theoretical Population Biology, 7, 256–276.
%Tajima, F. (1989) Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. Genetics, 123, 585–595.

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (nargin < 3)
    persite = 0;
end
if nargin < 2,
    fromS = 0; % estimate theta-W from S instead of m_mut.
end
if (isstruct(aln)), seq = aln.seq;
else seq = aln;
end

[n, m] = size(seq);
[S, ~, m_num] = countsegregatingsites(seq);
if (fromS)
    m_num = S;
end
nx = 1:(n - 1);
a1 = sum(1./nx);

t = m_num / a1;
if (nargout > 1),
    a2 = sum(1./nx.^2);
    tvar = t ./ a1 + (a2 * t^2) ./ a1^2;
end

if (persite),
    t = t / m;
end

%The classic "Watterson's Theta" statistic, generalized to missing data and multiple mutations per site:
%\[ \widehat\theta_w=\sum_{i=1}^{i=S}\frac{S}{\sum_{j=1}^{j=n_i-1}\frac{1}{j}} \]
%For this statistic, $S$ is either the number of segregating sites, or the number of mutations on the genealogy and $n_i$ is the sample size at site i. If totMuts == 1, the number of mutations is used, else the number ofsegregating sites is used.
% WTHETA - Watterson's theta
%
%Watterson, G. (1975) On the number of segragating sites in genetical models without recombination. Theoretical Population Biology, 7, 256–276.
%Tajima, F. (1989) Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. Genetics, 123, 585–595.

%    The classic "Watterson's Theta" statistic, generalized to missing data
%    and multiple mutations per site:
%    \f[
%    \widehat\theta_w=\sum_{i=1}^{i=S}\frac{S}{\sum_{j=1}^{j=n_i-1}\frac{1}{j}}
%    \f]\n
%    For this statistic, \f$S\f$ is either the number of segregating sites,
%    or the number of mutations on the genealogy and \f$n_i\f$ is the sample size
%    at site i. If totMuts == 1,
%    the number of mutations is used, else the number ofsegregating sites is used.
%    \warning Statistic undefined if there are untyped SNPs.  In the presence of
%    missing data, ThetaW is calculated as the sum (over all segregating sites)
%    of 1/a_sub_n, where a_sub_n is the denominator of ThetaW, using the number
%    of alleles without missing data as the sample size.  More formally, the
%    routine returns a calculation base on an unweighted sample size adjustment accross
%    sites.