function [status] = snp_writeplem(geno, filename)
%SNP_WRITEPLEM - saves as PLEM input file
%
% snp_writeplem(geno,filename)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (isempty(geno)), status = 0;
    return;
end
if (nargin < 2),
    [filename, pathname, filterindex] = uiputfile( ...
        {'*.txt', 'PLEM Input Format Files (*.txt)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Save as');
    if ~(filename), status = 0;
        return;
    end
    filename = [pathname, filename];
    if (filterindex == 1)
        if isempty(strfind(filename, '.'))
            filename = [filename, '.txt'];
        end
    end
end
fid = fopen(filename, 'wt');
if fid == -1
    status = 0;
    warning('SNP_WRITEPLEM:OpenFile', 'Unable to open file.');
    return;
end
g = snp_hhgeno(geno);
g(g == 3) = 0;
[n] = size(g, 1);

for k = 1:n
    fprintf(fid, '%d', g(k, :));
    fprintf(fid, '\n');
end
fclose(fid);

status = 1;
