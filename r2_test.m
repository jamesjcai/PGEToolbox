function [r2] = r2_test(aln,rvalue)
%R2_TEST - Ramos-Onsins & Rozas (R2) test
%
% Syntax: [r2] = r2test(aln,rvalue)
%
% Inputs:
%    aln   - Alignment structure or sequence matrix
%
% Outputs:
%    r2       - R_2 statistics
%    rvalue   - (optional) 3 for R3 test, 4 for R4 test
%
% Ramos-Onsins & Rozas (R2) test
% Detects recent severe population growth (indicated by low r2 value)
% This test is based on the difference between the number of singleton mutations and
% the average number of nucleotide differences. The R2 statistic is defined as
%
% \[
% R_2 = \frac{{({{\sum\limits_{i = 1}^n {(U_i  - k/2)} ^2 } \mathord{\left/
% {\vphantom {{\sum\limits_{i = 1}^n {(U_i  - k/2)} ^2 } n}} \right.
% \kern-\nulldelimiterspace} n})^{1/2} }}{S}
% \]
%
% where $n$ is the sample size, $S$ the total number of segregating sites, $k$ the
% average number of nucleotide differences between two sequences, and $U_i$ the
% number of singleton mutations in sequence $i$. The rationale of this test is that
% the expected numbers of singletons on a genealogy branch after a recent severe
% population growth event is $k/2$; consequently, lower values of $R_2$ are expected
% under this demographic scenario.
%
% Two R2 related tests namely, R3 and R4. These statistics differ from the R2 test
% in the power exponent values; in R3 and R4, the exponent values of 2 and 1/2 are
% replaced by 3 and 1/3, and by 4 and 1/4, respectively.
%
%REF: Ramos-Onsins SE, Rozas J.
%     Statistical properties of new neutrality tests against population growth.
%     Mol Biol Evol. 2002 Dec;19(12):2092-100.

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if (nargin<2), rvalue=2; end

if (isstruct(aln)), seq=aln.seq; else seq=aln; end
[n,m]=size(seq);
if (n<4) error('Four or more sequences are need to compute R_2 statistics'); end

[k] = thetapi(aln);
[S,V,m_num,sn_num,sm_num,Ui] = countsegregatingsites(aln);
%Ui is the vector of numbers of singleton mutations for each allele

r2=sum(abs(Ui-k/2).^rvalue);
r2=(r2/n)^(1/rvalue);
r2=r2/S;

if (nargout<1),
	i_dispheader('Ramos-Onsins and Rozas''s R2 Test')
	fprintf (['   R2 : %f\n'],r2);
	i_dispfooter
end
