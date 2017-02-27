function [Ds,Fs] = fuli93dsfs_simu(nsam,nrep,theta,segs,rho,nsites)
%FULI93DSFS_SIMU - Fu and Li's D* and F* statistics by coalescent simulation
%  Syntax: [Ds,Fs] = fuli93dsfs_simu(nsam,nrep,theta,segs,rho,nsites)

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
Ds=zeros(1,nrep);
Fs=zeros(1,nrep);

[G] = msrun(nsam,nrep,theta,segs,rho,nsites);
for k=1:nrep
        %io=ms(nsam,theta,segs);
	io=G{k};
	[segs]=size(io,2);
	x=sum(io,1);
	sigsegs=sum(x==1);
	xpi=sum(x.*(nsam-x));
	xpi=xpi/nn;
	[Ds(k),Fs(k)]=fuli93dsfs(nsam, segs, xpi, sigsegs);
end


if (nargout<1),
i_dispheader('Result of Fu and Li''s D*')
%	fprintf (['Average D* value: %f\n'],mean(Ds));
%	fprintf (['Average F* value: %f\n'],mean(Fs));
cireport(Ds)
i_dispheader('Result of Fu and Li''s F*')
cireport(Fs)
i_dispfooter
end
