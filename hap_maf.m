function [p, majalle, minalle] = hap_maf(haplo)
%HAP_MAF - Mininum Allele Frequence of markers in haplotype data
%[p] = hap_maf(haplodata)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

%n - individuls
%m - marker number

m = size(haplo, 2);
p = zeros(1, m);
majalle = 5 * ones(1, m);
minalle = 5 * ones(1, m);
for k = 1:m
    x = haplo(:, k);
    x(x == 5) = [];
    if ~isempty(x)
        [a, ~, c] = unique(x);
        if length(a) == 1
            p(k) = 0;
            majalle(k) = a;
            minalle(k) = a;
        else
            y = sum(c == 1) / length(x);
            [y, idx] = min([1 - y, y]);
            p(k) = y;
            majalle(k) = a(idx);
            b = setdiff(a, a(idx));
            minalle(k) = b(1);
        end
    else
        p(k) = nan;
        majalle(k) = nan;
        minalle(k) = nan;
    end
end
