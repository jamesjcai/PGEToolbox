function [thew] = thetaw_simu(nsam,nrep,theta,segs,rho,nsites)
%THETAW_SIMU - Watterson's theta by coalescent simulation
% Syntax: [thew] = thetaw_simu(nsam,nrep,theta,segs,rho,nsites)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (nargin<6), nsites=2; end
if (nargin<5), rho=0; end
if (nargin<4), segs=0; end
if (nargin<2), nrep=1000; end
if (nargin<1), error('Need sample size.'); end

nn=nchoosek(nsam,2);
thew=zeros(1,nrep);


[G] = msrun(nsam,nrep,theta,segs,rho,nsites);
a1=sum(1./(1:(nsam-1)));
for (k=1:nrep),
      io=G{k};
      thew(k)=size(io,2)/a1;
end

if (nargout<1),
	i_dispheader('Watterson''s theta (theta_W)')
	cireport(thew)
	i_dispfooter
end

