function [h] = faywu00h_test(aln,ancseq,fromS)
%FAYWU00H_TEST - Fay and Wu's test
%  Syntax: [h] = faywu00h_test(aln,ancseq,fromS)
%
% Test statistic, H, as the difference between $\theta_pi$, which is
% influenced most by variants at intermediate frequencies, and $\theta_H$,
% which is influenced most by high-frequency variants.

% This test is analogous to and was motivated by the F(0, -1) test described in FU
% 1995 , the only difference being F(0, -1) is the difference between {pi} and H
% divided by the variance of the difference. Under neutrality the expected
% difference between two estimators of {theta} is zero, but following a
% hitchhiking event theta_H and theta_W should be >theta_{pi}.

% ... a statistic, H, to measure an excess of high compared to
% intermediate frequency variants. Only a few high-frequency variants are
% needed to detect hitchhiking since not many are expected under neutrality

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<3,
    fromS=0;    % estimate theta-W from S instead of m_mut.
end
if nargin<2,
    ancseq=0;
end

if (isstruct(aln)), seq=aln.seq; else seq=aln; end
m = size(seq,2);
if (ancseq)
	[thissfs,smpsize]=sfs(seq,ancseq);
else
	[thissfs,smpsize]=sfs(seq);
end

%smpsize=n-(max(size(ancseq))==1);

[S,V,m_num,sn_num,sm_num] = countsegregatingsites(aln);

if (fromS), Sk=S; else Sk=m_num; end
h=[]; hraw=[];
if ~(Sk==0&&sum(thissfs(:))==0)
    [h,hraw] = faywu00h(smpsize, Sk, thissfs);
end

%thetapi(aln,0,1)

if (nargout<1),
i_dispheader('Fay and Wu''s H Test')
	fprintf('No. of Sequences (n) : %d\n', smpsize);
	fprintf('No. of Sites (m) : %d\n', m);
	fprintf('\n');
	%fprintf(['Recombination rate (c) : %f\n'], rcbrate);
	%fprintf(['Probability of back-mutation : %f\n'], bckrate);
	%fprintf(['\n']);
    
    if isempty(h)
        fprintf ('Fay and Wu''s H = -NULL-\n');
        fprintf ('Fay and Wu''s H (normalized) = -NULL-\n');    
    else
        fprintf ('Fay and Wu''s H = %f\n', hraw);
        fprintf ('Fay and Wu''s H (normalized) = %f\n', h);
        [thew] = thetaw(aln);
        [H] = faywu00h_simu(smpsize,10000,thew,0);
        p=2.*sum(H>abs(hraw))./10000;
        %p=sum(hraw<H)./1000;
        fprintf ('Statistical significance:\n P = %f%s\n',p,sigtag(p));
    end
i_dispfooter
end




% The input parameters are:
%
% Sample size = 10
%
% Recombination rate = 0
%
% Divergence = 0
%
% Loops = 5000
%
% The results from the htest are given below:
%
% Observed data: 3 2 1 0 0 0 0 1 1
% Back-mutation: 0.000000
% Observed: Theta(w) = 2.827886, Theta(pi) = 2.333333, Theta(h) = 3.666667
% Probability: D = 0.256400 H = 0.154200