function reportpolysites(aln)
%REPORTPOLYSITES - Report polymorphic sites
%An example report:
% Number of sequences: 4    Number of sequences used: 4
% Selected region: 1-55   Number of sites: 55
%    Number of sites (excluding fixed gaps / missing data): 48
% Total number of sites (excluding sites with gaps / missing data): 46
%   Sites with alignment gaps or missing data: 9
%   Invariable (monomorphic) sites: 36
%   Variable (polymorphic) sites: 10   (Total number of mutations: 11)
%      Singleton variable sites: 9
%      Parsimony informative sites: 1
%
%      Singleton variable sites (two variants): 8
%         Site positions:  11 17 21 31 33 37 40 44
%      Parsimony informative sites (two variants): 1
%         Site positions:  27
%      Singleton variable sites (three variants): 1
%         Site positions:  20
%      Parsimony informative sites (three variants): 0
%      Variable sites (four variants): 0
%
% Protein Coding Region assignation:  No
%
% Syntax: reportpolysites(aln)
%
% Inputs:
%    aln   - Alignment structure
%
% See also:

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


[N,M]=size(aln.seq);

Aln2=rmgaps(aln);
[n,m]=size(Aln2.seq);
if (n<4) error('n must > 4'); end

	i_dispheader('Statistics of Polymorphic Sites')
	disp('Mode : Nucleotides');
	fprintf(['No. of Sequences (n) : %d\n'],n);
	fprintf(['No. of Sites (m) : %d\n'],M);
	fprintf(['Total number of sites\n (excluding sites with gaps/missing data) : %d\n'],m);
	fprintf(['\n']);
	fprintf(['Sites with alignment gaps or missing data : %d\n'],M-m);
	[Mono, v_site]=countinvariablesites(aln);

	fprintf(['Invariable (monomorphic) sites : %d\n'],Mono);

	[s_site,v_site,m_num,s_num]=countsegregatingsites(aln);

	fprintf(['Variable (polymorphic) sites : %d\n'],s_site);
	fprintf(['   Total number of mutations: %d\n'],m_num);
	fprintf(['   Singleton variable sites : %d\n'], s_num);

	[Aln3,sites]=extractinformativesites(aln);
	fprintf(['   Parsimony informative sites : %d\n\n'], size(sites,2));
	i_dispfooter