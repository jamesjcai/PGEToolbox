function [status] = snp_writegenepop(genodata, filename, subpopidx)
%snp_writegenepop - save genotype data for GenePop
% [status] = snp_writegenepop(genodata,filename,subpopidx)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

status = 0;
[smpln, markn] = snp_samplen(genodata);
smpln = smpln / 2;

if nargin < 3
    subpopidx = ones(1, smpln);
    names = {};
    for (k = 1:smpln),
        names{k} = ['Idv_', num2str(k)];
    end
    [s, v] = choosebox('Name', 'Separate individuals', 'PromptString', ...
        'Subpopulation 1:', 'SelectString', 'Subpopulation 2:', ...
        'ListString', names);
    subpopidx(s) = 2;
    if ~(v == 1), return; end
end


if nargin < 2
    [filename, pathname, filterindex] = uiputfile( ...
        {'*.*'}, 'Save as');
    if ~(filename), status = 0;
        return;
    end
    filename = [pathname, filename];
end

if strfind(filename, '.') > 0
    fid1 = fopen(filename, 'w');
else
    fid1 = fopen([filename, '.dat'], 'w');
end

%pos=markinfo.pos;
%rsid=markinfo.rsid;

%First line: any character. Use this line to store information about your
%data.
fprintf(fid1, 'GenePop input file (%d pop, %d loci)\n', max(subpopidx(:)), markn);
for k = 1:markn
    fprintf(fid1, 'loc_%d\n', k);
end


oldpop = 0;
indc = 1;
for i = 1:smpln
    if subpopidx(i) ~= oldpop
        fprintf(fid1, 'pop\n');
        oldpop = subpopidx(i);
    end
    fprintf(fid1, 'ind%d, ', indc);
    indc = indc + 1;
    for j = 1:2:markn * 2
        if (genodata(i, j) == 5 || genodata(i, j+1) == 5)
            fprintf(fid1, ' 0000 ');
        else
            fprintf(fid1, '0%d0%d ', genodata(i, j), genodata(i, j+1));
        end
    end
    fprintf(fid1, '\n');
end

fclose(fid1);
status = 1;
