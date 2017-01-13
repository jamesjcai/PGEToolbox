function [h,hvar] = thetah(aln)
%THETAH - Fay's theta_H based on the frequency of derived alleles
%
% Syntax: [theta_h] = thetah(aln)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (isstruct(aln)), seq=aln.seq; else seq=aln; end
[freqall,~] = sfs(seq);
n=size(seq,1);

%warning('The first sequence in alignment was used as ancestor sequence. Number of sequences becomes n-1.');
nx=1:(n-1);
h=2*sum((nx.^2).*freqall)./(n*(n-1));

%old error
%a1 = sum(1./nx);
%theta_h=sum(2*freqall.*nx.^2)/a1;

%[theta_pi] = thetapi(seq,0,0)     % not per site, but using sfs

if nargout>1
an=sum(1./nx);
bn=sum(1./(nx.^2));
bn2=sum(1./((1:n).^2));

%[S,V,m_num,s_num,sm_num] = countsegregatingsites(seq);
S=sum(freqall);

%In practice, {theta} in (12) can be estimated by {theta}W, and
%{theta}2 can be estimated by s(s-1)/(an^2+bn) (TAJIMA 1989).
t1=thetaw(seq);
t2=S*(S-1)./(an.^2+bn);
%Equation 10
%http://www.genetics.org/cgi/content/full/174/3/1431

hvar=t1+2*(36*(n^2)*(2*n+1)*bn2-116*n^3+9*n^2+2*n-3)*t2/(9*n*(n-1)^2);
end
