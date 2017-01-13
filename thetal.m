function [thetal,thetalvar] = thetah(aln)
%THETAL - Zeng's theta_L based on the frequency of derived alleles
%
% Syntax: [l,lvar] = thetal(aln)
%
%REF:
%Zeng K, Fu YX, Shi S, Wu CI.
%Statistical tests for detecting positive selection by utilizing high-frequency variants.
%Genetics. 2006 Nov;174(3):1431-9. PMID: 16951063

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if (isstruct(aln)), seq=aln.seq; else seq=aln; end
[freqall,nsampl] = sfs(seq);
[n,m]=size(seq);
nx=1:(n-1);
thetal=sum(nx.*freqall)./(n-1);


%[theta_pi] = thetapi(seq,0,1);     % not per site, but using sfs

if nargout>1
an=sum(1./nx);
bn=sum(1./(nx.^2));
bn2=sum(1./([1:n].^2));

%[S,V,m_num,s_num,sm_num] = countsegregatingsites(seq);
S=sum(freqall);

%In practice, {theta} in (12) can be estimated by {theta}W, and
%{theta}2 can be estimated by s(s-1)/(an^2+bn) (TAJIMA 1989).
t1=thetaw(seq);
t2=S*(S-1)./(an.^2+bn);
thetalvar=t1*n/(2*(n-1))+((2*(n/(n-1))^2)*(bn2-1)-1)*t2;
end
