function [f] = watterson78f(aln)
%WATTERSON78F - Watterson's homozygosity test of neutrality (1978)
% Syntax: [f] = watterson78f(aln)
%
% The homozygosity statistic is the sum of squared allele frequencies
% (proportions). The statistic can be used to test deviation from
% neutrality expectations.
%
%REF:  Watterson, G. A.
%THE HOMOZYGOSITY TEST OF NEUTRALITY
%Genetics 1978 88: 405-417
%http://www.genetics.org/cgi/reprint/88/2/405.pdf
%
% Ewens, W.J. 1972. The sampling theory of selectively neutral alleles.
% Theor. Popul. Biol. 3: 87–112.
%
% See also: HAPDIV_TEST,

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


% haplotype homozygosity conditional on
% the number of haplotypes) as the Ewens-Watterson (EW) test
% statistic (Watterson 1978).

if (isstruct(aln)), seq=aln.seq; else seq=aln; end
%[n,m]=size(seq);

[~,sizHap]=counthaplotype(seq);
p=sizHap./sum(sizHap);
f=sum(p.^2);    % Equation on page 411

if (nargout<1),
i_dispheader('Ewens-Watterson''s Homozygosity Test')
	fprintf('Haplotype homozygosity, F = %f\n',f);
i_dispfooter
end

% P-value < 0.025: Too even -> Balancing selection or recent bottleneck
% P-value > 0.975: Too uneven -> Directional selection or population growth

%Watterson’s test is conditioned on the number of haplotypes, whereas Depaulis and Veuille's H
%is conditioned on S (and n).



% A subsequent “exact test” was based on whether
% C itself was unlikely given n and K (Slatkin 1994, 1996).

% REF: http://www.genetics.org/cgi/content/short/169/3/1763

