function [status] = snp_writesnphap(geno, filename)
%SNP_WRITESNPHAP - saves as SNPHAP input file format
%
% snp_writesnphap(geno)
% snp_writesnphap(geno,filename)

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
        {'*.snphap;*.dat', 'SNPHAP Format Files (*.snphap, *.dat)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Save as');
    if ~(filename), status = 0;
        return;
    end
    filename = [pathname, filename];
    if (filterindex == 1)
        if (isempty(find(filename == '.')))
            filename = [filename, '.dat'];
        end
    end
end
fid = fopen(filename, 'wt');
if (fid == -1),
    status = 0;
    warning('Unable to open file.');
    return;
end

[n, m] = size(geno);
%[geno] = i_simpgeno(geno);
[geno] = snp_12geno(geno);

for (k = 1:n),
    fprintf(fid, ['%d\t'], k);
    for (j = 1:m - 1),
        fprintf(fid, ['%d\t'], geno(k, j));
    end
    fprintf(fid, ['%d\n'], geno(k, j+1));
end
fclose(fid);
status = 1;


    function [geno2] = i_simpgeno(geno)
        %SIMPGENO - Simplify GENO coding convention

        geno2 = [];
        n = snp_marklen(geno);
        for k = 1:n
            gen = geno(:, 2*k-1:2*k);
            genx = gen(gen < 5);
            [a] = unique(genx);
            if length(a) == 1
                gen(:) = 1;
            elseif length(a) == 2
                if (sum(genx == a(1)) >= sum(genx == a(2)))
                    gen(gen == a(1)) = 1;
                    gen(gen == a(2)) = 2;
                else
                    gen(gen == a(1)) = 2;
                    gen(gen == a(2)) = 1;
                end

            end
            gen(gen > 4) = 0;
            geno2 = [geno2, gen];
        end
