function [D] = hapdiv_simu(nsam, nrep, theta, segs, rho, nsites)
%HAPDIV_SIMU - Haplotype diversity/heterozygosity by coalescent simulation

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
D = zeros(1, nrep);

[G] = msrun(nsam, nrep, theta, segs, rho, nsites);
for (k = 1:nrep),
    %io=ms(nsam,theta,segs);
    io = G{k};
    [numHap, sizHap] = counthaplotype(io);
    D(k) = hapdiv(sizHap, 1);
end


if (nargout < 1),
    i_dispheader('Haplotype diversity/heterozygosity D')
    %	fprintf (['Average value: %f\n'],mean(D));
    cireport(D)
    i_dispfooter
end
