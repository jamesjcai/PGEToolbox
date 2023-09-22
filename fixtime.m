function [t] = fixtime(p, s, Ne)
%FIXTIME - time to fixation (generations)
%
%Usage: [u]=fixtime(p,s,Ne)
%
% Kimura & Ohta 1969

% N=1000;
% s=linspace(-0.04,0.04,50)
% p=1/(2*N);

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin < 3, Ne = 1000; end
if nargin < 2, s = 0.01 ./ Ne; end
if nargin < 1, p = 1 ./ Ne; end

S = 2 * Ne .* s;

c = 2; % c=1 for haploids, c=2 for diploids

%F = @(x) (exp(-2*S.*x)-1).*(exp(-2*S.*(1-x))-1)./(s.*x.*(1-x).*(exp(2*S)-1));
F = @(x) (1 - exp(-2*S.*x)) .* (1 - exp(-2*S.*(1 - x))) ./ (s .* x .* (1 - x) .* (1 - exp(2*S)));

t = i_quadl(F, p);


    function y = i_quadl(F, x)
        n = length(x);
        y = zeros(1, n);
        for k = 1:n
            y(k) = quadl(F, x(k), 1);
        end
