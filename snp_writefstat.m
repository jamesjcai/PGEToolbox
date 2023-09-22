function [status] = snp_writefstat(genodata, filename, subpopidx)
%snp_writefstat - save genotype data for FSTAT
% [status] = snp_writefstat(genodata,filename,subpopidx)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

status = 0;

[indvn, m] = size(genodata);
markn = m ./ 2;

if (markn > 200)
    warning on
    warning('Fstat may not allow markers > 200')
    warning off
end

if nargin < 3
    subpopidx = ones(1, indvn);
    names = {};
    for (k = 1:indvn),
        names{k} = ['Idv_', num2str(k)];
    end
    [s, v] = choosebox('Name', 'Separate individuals', 'PromptString', ...
        'Subpopulation 1:', 'SelectString', 'Subpopulation 2:', ...
        'ListString', names);
    subpopidx(s) = 2;
    if (v ~= 1), return; end
end


if nargin < 3
    [filename, pathname, filterindex] = uiputfile( ...
        {'*.*'}, 'Save as');
    if ~(filename), status = 0;
        return;
    end
    filename = [pathname, filename];
end

fid1 = fopen([filename, '.dat'], 'w');

%pos=markinfo.pos;
%rsid=markinfo.rsid;


fprintf(fid1, '   2  %d  4  1\n', markn);
for k = 1:markn
    fprintf(fid1, 'loc_%d\n', k);
end

for i = 1:indvn
    fprintf(fid1, '%d   ', subpopidx(i));
    for j = 1:2:markn * 2
        if (genodata(i, j) == 5 || genodata(i, j+1) == 5)
            fprintf(fid1, ' 0 ');
        else
            fprintf(fid1, '%d%d ', genodata(i, j), genodata(i, j+1));
        end
    end

    fprintf(fid1, '\n');
end

fclose(fid1);
status = 1;
