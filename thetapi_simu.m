function [thepi] = thetapi_simu(nsam,nrep,theta,segs,rho,nsites)
%THETAH_SIMU - Fay's theta_H statistic by coalescent simulation
% Syntax: [thepi] = thetapi_simu(nsam,nrep,theta,segs,rho,nsites)

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
thepi=zeros(1,nrep);

[G] = msrun(nsam,nrep,theta,segs,rho,nsites);
for (k=1:nrep),
        io=G{k};
	x=sum(io,1);
	thepi(k)=2*sum(x.*(nsam-x))/(nsam*(nsam-1));
    %thepi(k)=sum(x.*(nsam-x))./nchoosek(nsam,2);
end

if (nargout<1),
	i_dispheader('Result of Nei''s theta (pi)')
	% fprintf (['Average value: %f\n'],mean(H));
	cireport(thepi)
	i_dispfooter
end

