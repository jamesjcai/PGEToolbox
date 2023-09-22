function [H] = faywu00h_simu(nsam, nrep, theta, segs, rho, nsites)
%FAYWU00H_SIMU - Fay and Wu's H statistic by coalescent simulation
%  Syntax: [H] = faywu00h_simu(nsam,nrep,theta,segs,rho,nsites)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (nargin < 6), nsites = 2; end
if (nargin < 5), rho = 0; end
if (nargin < 4), segs = 0; end
if (nargin < 2), nrep = 1000; end
if (nargin < 1), error('Need sample size.'); end

nn = nchoosek(nsam, 2);
H = zeros(1, nrep);

[G] = msrun(nsam, nrep, theta, segs, rho, nsites);
for k = 1:nrep
    io = G{k};
    [segs] = size(io, 2);

    x = sum(io, 1);

    %%xpi=sum(x.*(nsam-x));
    %%xh=sum(x.^2);
    %%xfayh=xpi-xh;
    %xfayh=sum(x*nsam-x.^2);  % this is a simplified form of above three
    %H(k)=xfayh/nn;

    %/* Fay's theta_H  */
    %thetah=2*sum(x.^2)/(nsam*(nsam-1));
    p2 = 2 * (x ./ nsam);
    H(k) = -sum(p2.*(p2 - 1)*(nsam / (nsam - 1)));

end

if nargout < 1
    i_dispheader('Result of Fay and Wu''s H')
    %fprintf (['Average value: %f\n'],mean(H));
    cireport(H)
    i_dispfooter
end
