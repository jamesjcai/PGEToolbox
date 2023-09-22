function [f, f1] = stnpdfsmpl(n, h, s, Ne)
%
%[f1,f1a]=stnpdfsmpl(10,.5,-2)
%[f2,f2a]=stnpdfsmpl(10,.5,-5)
%bar([f1;f2]')
%legend({'2N_{e}s=-2','2N_{e}s=-5'})
%
% equ (2)
% http://www.plosgenetics.org/article/info:doi/10.1371/journal.pgen.0030163
%
%components of the site-frequency spectrum (X1, X2, . . ., Xn) as
%independent Poisson random variables with mean:
%E(Xi)=theta/2*int(binopdf(i,n,x).*stnpdf(x,Sx);)
%
% See also: stnpdf, stnpdfh, stnpdfsmpl_test

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin == 4, Sx = 2 * Ne * s; end
if nargin == 3, Sx = s; end
if isempty(h), h = 0.5; end
semidom = 0;
if (h == 0.5), semidom = 1; end % semidom

f1 = zeros(1, n-1);
for i = 1:n - 1;
    if (semidom) % semidom
        F = @(x) binopdf(i, n, x) .* stnpdf(x, Sx);
    else
        F = @(x) binopdf(i, n, x) .* stnpdfh(x, h, Sx);
    end
    f1(i) = quadl(F, 0, 1);
end

f = f1 ./ sum(f1);

%2*(1-exp(-Sx.*(1-p)))./((1-exp(-Sx)).*(p.*(1-p)))
