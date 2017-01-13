function [theta] = estimatetheta(aln,id)
%ESTIMATETHETA - Estimate theta

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (nargin<2), id=1; end
if (isstruct(aln)), seq=aln.seq; else seq=aln; end

if (nargout<1),
    [N,M]=size(seq);
	[theta1,theta1var]=thetapi(seq);
	[theta2,theta2var]=thetaw(seq);
	[theta3,theta3var]=thetaw(seq,1);
 %   [theta4,theta4var]=thetah(seq);

	%[theta3]=thetaEwen(seq);
	i_dispheader('Estimation of theta (i.e., 4Nu)')
	fprintf (['Theta_W (from Eta*) = %f, var = %f (per sequence)\n'], theta2,theta2var);
	fprintf (['                    = %f (per site)\n'], theta2/M);
	fprintf (['Theta_W (from S**)  = %f, var = %f (per sequence)\n'], theta3,theta3var);
	fprintf (['                    = %f (per site)\n'], theta3/M);
	fprintf (['Theta_Pi***         = %f, var = %f (per sequence)\n'], theta1,theta1var);
	fprintf (['                    = %f (per site)\n'], theta1/M);
%	fprintf (['Theta_H***         = %f, var = %f (per sequence)\n'], theta4,theta4var);
%	fprintf (['                    = %f (per site)\n'], theta4/M);
%   fprintf (['Ewens''s estimate (formula 9.26) = %f\n'], theta3);
	fprintf ('\n');
	fprintf ('*Eta is the total number of mutations\n');
	fprintf ('**S is number of segregating sites\n');
	fprintf ('***Theta_Pi equals nucleotide diversity\n');

    i_dispfooter

else
	switch (id)
	    case (1)
		 theta=thetapi(seq);
	    case (2)
		 theta=thetaw(seq);
	    case (3)
		 theta=thetaewens(seq);
	end
end

%thetah(aln);
%Calculate Theta ( = 4Nu) from site homozygosity, a la Fay and Wu (2000). This
%statistic is problematic in general to calculate when there are multiple hits.
%The test requires that the ancestral state (inferred from the outgroup) still
%be segregating in the ingroup. If that is not true, the site is skipped.

%thetal(aln);
%Calculate Theta ( = 4Nu) from site homozygosity, corresponding to equation 1 in
%Thornton and Andolfatto (Genetics) "Approximate Bayesian Inference reveals evidence
%for a recent, severe, bottleneck in a Netherlands population of Drosophila
%melanogaster," (although we labelled in $theta_\eta$ in that paper) The test
%requires that the ancestral state (inferred from the outgroup) still be segregating
%in the ingroup. If that is not true, the site is skipped.
