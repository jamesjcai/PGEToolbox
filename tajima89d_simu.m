function [D] = tajima89d_simu(nsam, nrep, theta, segs, rho, nsites)
%TAJIMA89D_SIMU - Tajima's D statistic by coalescent simulation
% Syntax: [D] = tajima89d_simu(nsam,nrep,theta,segs,rho,nsites)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2014-03-26 09:30:26 -0500 (Wed, 26 Mar 2014) $
% $LastChangedRevision: 758 $
% $LastChangedBy: jcai $

if (nargin < 6), nsites = 2; end
if (nargin < 5), rho = 0; end
if (nargin < 4), segs = 0; end
if (nargin < 2), nrep = 1000; end
if (nargin < 1), error('Need sample size.'); end

nn = nchoosek(nsam, 2);
D = zeros(1, nrep);

[G] = msrun(nsam, nrep, theta, segs, rho, nsites);
for k = 1:nrep
    %io=ms_mex(nsam,theta,segs);
    io = G{k};
    if isempty(io)
        D(k) = 0;
    else
        [segs] = size(io, 2);

        x = sum(io, 1);
        xpi = sum(x.*(nsam - x));
        %xh=sum(x.^2); xfayh=sum(x*nsam-x.^2);
        xpi = xpi / nn;
        %xh=xh/nn; xfayh=xfayh/nn;

        % xpi = sum(sum(io*(1-io)'+(1-io)*io'))/nn/2;  Daniel Powell's
        % homework answer

        D(k) = tajima89d(nsam, segs, xpi);
    end
end

if (nargout < 1),
    i_dispheader('Result of Tajima''s D')
    %	fprintf (['Average value: %f\n'],mean(D));
    cireport(D)
    i_dispfooter
end
