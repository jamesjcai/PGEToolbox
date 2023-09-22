function [ihs] = snp_ihsregion(chrid, startn, endn, radius, popcode)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin < 5
    popcode = 'CEU';
else
    popcode = upper(popcode);
    if ~(ismember(popcode, {'CEU', 'CHB', 'JPT', 'YRI', 'JPT+CHB'})),
        error('Not available population code.');
    end
    if (strcmp(popcode, 'JPT+CHB'))
        popcode = 'JPT%2BCHB';
    end
end

if nargin < 4
    radius = 1000000;
end


query = sprintf('Chr%d:%d..%d', chrid, startn, endn);
disp(sprintf('Dowloading phased haplotype data from target region %s...', query))
[s0, hapldata1, markinfo1] = snp_downloadhaplotype(query, popcode);

if isempty(hapldata1)
    ihs = [];
    return;
end


startn2 = max(1, startn-radius);
endn2 = min(chrlen(chrid), endn+radius);


query = sprintf('Chr%d:%d..%d', chrid, startn2, endn2);
disp(sprintf('Dowloading phased haplotype data from extend region %s...', query))
[s0, hapldata2, markinfo2] = snp_downloadhaplotype(query, popcode);

[yes, idxx] = ismember(markinfo1.rsid, markinfo2.rsid);
ihs = zeros(1, length(idxx));

for k = 1:length(idxx)
    idx = idxx(k);
    disp(sprintf('%d of %d ...', k, length(idxx)))
    ehh = snp_ehh(hapldata2, idx);
    if (size(ehh, 1) == 2),
        ihh1 = trapz(ehh(1, :));
        ihh2 = trapz(ehh(2, :));
        if ihh2 > 0
            ihs(k) = abs(log(ihh1/ihh2));
        end
    end
end
