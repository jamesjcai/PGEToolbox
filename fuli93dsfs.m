function [ds,fs] = fuli93dsfs(nsam,segs,neipi,sigsegs)
%FULI93DSFS - Compute Fu and Li (1993)'s D* and F* statistics
%
% Syntax: [ds,fs] = fuli93dsfs(nsam,segs,neipi,sigsegs)
%
% Inputs:
%    nsam      - number of samples
%    segs      - number of segregating sites
%    neipi     -
%    sigsegs   -
%
% Outputs:
%    ds   - D*
%    fs   - F*
%
% See also:

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


n=nsam; m_num=segs; k=neipi; sm_num=sigsegs;

nx=1:(n-1);
a1 = sum(1./nx);
a2 = sum(1./nx.^2);

% The following coefficients are for Fu and Li d* and f*
if (n==2),
	cn=1;
else
	cn=2*(((n*a1)-(2*(n-1)))/((n-1)*(n-2)));
end   %% FU LI eq 14. p 695

anp1=a1+(1/n);
dn=(cn+((n-2)/((n-1)^2))+((2/(n-1))*(3/2-(((2*anp1)-3)/(n-2))-(1/n))));
v_Dst=((((n/(n-1))^2)*a2)+((a1^2)*dn)-((2*n*a1*(a1+1))/((n-1)^2)))/((a1^2)+a2);
u_Dst=((n/(n-1))*(a1-(n/(n-1))))-v_Dst;

% Stat is Fu and Li D*
stat_numDs = (m_num * n / (n - 1)) - (sm_num * a1);
stat_denDs = sqrt ( (u_Dst * m_num) + (v_Dst * (m_num ^2)));
stat_denDs=max(stat_denDs,eps);
ds=stat_numDs/stat_denDs;

if (nargout>1),
	% The below (corrected) expression is from Simonsen et al p 428
	v_Fst=(((2*(n^3)+(110*(n^2))-(255*n)+153)/(9*(n^2)*(n-1)))...
			+((2*(n-1)*a1)/(n^2))-(8*a2/n))/((a1^2)+a2);
	u_Fst=((((4*(n^2))+(19*n)+3-(12*(n+1)*anp1))/(3*n*(n-1)))/a1)-v_Fst;

	% Stat is Fu and Li F*
	stat_numFs = k - (sm_num * (n - 1) / n);
	stat_denFs = sqrt ( (u_Fst * m_num) + (v_Fst * (m_num ^ 2)));
	stat_denFs=max(stat_denFs,eps);
	fs=stat_numFs/stat_denFs;
end
