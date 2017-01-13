function [h] = zeng06e_test(aln,fromS)
% Zeng's E Test
% [h] = zeng06e_test(aln,fromS)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<2,
    fromS=0;    % estimate theta-W from S instead of m_mut.
end

if (isstruct(aln)), seq=aln.seq; else seq=aln; end
[n,m] = size(seq);
[thissfs,smpsize]=sfs(seq);

%smpsize=n-(max(size(ancseq))==1);

[S,V,m_num,sn_num,sm_num] = countsegregatingsites(aln);


if (fromS),
    [h,hraw] = zeng06e(smpsize, S, thissfs);
else
    [h,hraw] = zeng06e(smpsize, m_num, thissfs);
end

%thetapi(aln,0,1)

if (nargout<1),
i_dispheader('Zeng et al''s E Test')
	fprintf(['No. of Sequences (n) : %d\n'], smpsize);
	fprintf(['No. of Sites (m) : %d\n'], m);
	fprintf(['\n']);
	%fprintf(['Recombination rate (c) : %f\n'], rcbrate);
	%fprintf(['Probability of back-mutation : %f\n'], bckrate);
	%fprintf(['\n']);
	fprintf (['Zeng et al''s E = %f\n'], hraw);
	fprintf (['Zeng et al''s H (normalized) = %f\n'], h);
i_dispfooter;
end
