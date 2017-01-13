function [Ds,Fs] = fuli93dsfs_test(aln,switchoutput)
%FULI93DSFS_TEST - Fu and Li's D* and F* tests of neutrality
% D* test statistic
% The D* test statistic is based on the differences between hs, the number
% of singletons (mutations appearing only once among the sequences), and h,
% the total number of mutations (Fu and Li 1993, p. 700 bottom).
%
% F* test statistic
% The F* test statistic is based on the differences between hs, the number of
% singletons (mutations appearing only once among the sequences), and k, the
% average number of nucleotide differences between pairs of sequences (Fu and Li
% 1993, p. 702; see also Simonsen et al. 1995, equation 10).
%
% Syntax: [Ds,Fs]=fuli93dsfs_test(aln)
%
% Inputs:
%    aln   - Alignment structure
%
% Outputs:
%    Ds     - Fu and Li's D* statistic
%    Fs     - Fu and Li's F* statistic
%
%
% See also:

% Fu Y.X and Li W.H. (1993) "Statistical Tests of Neutrality of
% Mutations." Genetics 133:693-709.
%
% Fu Y.X. (1996) "New Statistical Tests of Neutrality for DNA samples
% from a Population." Genetics 143:557-570.
%
% Tajima F. (1989) "Statistical method for testing the neutral mutation
% hypothesis by DNA polymorphism." Genetics 123:585-595.


% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (nargin<2), switchoutput=0; end

if (isstruct(aln)), seq=aln.seq; else seq=aln; end

[n,m]=size(seq);
if n<4
    error('Four or more sequences are need to compute and Fu and Li''s statistics'); 
end

[k] = thetapi(aln);
[S,~,m_num,~,sm_num] = countsegregatingsites(aln);



[Ds,Fs]=fuli93dsfs(n,m_num,k,sm_num);
%[Ds,Fs]=fuli93dsfs(n,S,k,sn_num);

% double SequenceStatistics::fuliD(const PolymorphismSequenceContainer & ingroup, const PolymorphismSequenceContainer & outgroup) {
% 	unsigned int n = ingroup.getNumberOfSequences();
% 	double nn = (double) n;
% 	map<string, double> values = _getUsefullValues(n);
% 	double vD = 1. + (pow(values["a1"], 2) / (values["a2"] + pow(values["a1"], 2))) * (values["cn"] - ((nn + 1.) / (nn - 1.)));
% 	double uD = values["a1"] - 1. - vD;
% 	double eta = (double) totNumberMutations(ingroup);//using the number of mutations
%         //double eta = (double)polymorphicSiteNumber(ingroup);
% 	double etae = (double) totMutationsExternalBranchs(ingroup,outgroup);
% 	return (eta - values["a1"] * etae) / sqrt((uD * eta) + (vD * eta * eta));
% }

% % Stat is Fu and Li D
% a=a1; b=a2;
% c=2*(n*a-2*(n-1))/((n-1)*(n-2));
% v=1+a*a/(b+a*a)*(c+(n+1)/(1-n));
% u=a-1.0-v;
% % // etae = totMutationsExternalBranchs(ingroup,outgroup);
% stat_numD = m_num-a*etae;
% stat_denD = sqrt ( (u * m_num) + (v * (m_num ^2)));
% D=stat_numD/stat_denD;
%
% % Stat is Fu and Li F
% stat_numF = k - sm_num;
% stat_denF = sqrt ( (u_Fst * m_num) + (v_Fst * (m_num ^ 2)));
% F=stat_numF/stat_denF;



if (nargout<1),
i_dispheader('Fu and Li''s D* and F* Neutrality Test')
	disp('Mode : Nucleotides');
	disp('Gaps/Missing data : Complete Deletion');

	fprintf('No. of Sequences (n) : %d\n', n);
	fprintf('No. of Sites (m) : %d\n', m);
	fprintf('\n');
	fprintf('No. of Segregating sites (S) : %d\n', S);
	fprintf('Total number of mutations, (Eta) : %d\n', m_num);
	fprintf('\n');
	fprintf('Average number of nucleotide differences : %f\n', k); % scaled mutation rate theta (k)
	fprintf('Nucleotide diversity (per site), Pi : %f\n', nucdiv(aln));
	fprintf('\n');
	%fprintf (['Diff = %f, s.e. = %f\n'], Diff, DiffSE);

	fprintf ('Fu and Li''s D* = %f\n', Ds);

	[theta] = thetaw(aln);
	[Ds_rep,Fs_rep] = fuli93dsfs_simu(n,1000,theta,0);

    p=2.*sum(Ds_rep>abs(Ds))./1000;
	fprintf ('Statistical significance:\n P = %f%s\n',p,sigtag(p));
	p=max(sum(Fs>Fs_rep),sum(Fs<Fs_rep))./1000;
	fprintf ('Fu and Li''s F* = %f\n', Fs);
	fprintf ('Statistical significance:\n P = %f%s\n',p,sigtag(p));
i_dispfooter
end

if (switchoutput),
    temp=Ds; Ds=Fs; Fs=temp;
end