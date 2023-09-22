function [hvar] = thetavar(n, S, typeid)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

nx = 1:n - 1;
an = sum(1./nx);
bn = sum(1./(nx.^2));
bn2 = sum(1./([1:n].^2));

t1 = S ./ an;
t2 = S .* (S - 1) ./ (an.^2 + bn);


switch typeid

    case 1 %thetaH
        hvar = t1 + 2 .* (36 .* (n.^2) .* (2 .* n + 1) .* bn2 - 116 .* n.^3 + 9 .* n.^2 + 2 .* n - 3) .* t2 ./ (9 .* n .* (n - 1).^2);
    case 2 %thetaW
        hvar = t1 ./ an + (bn * t1.^2) ./ an.^2;
end