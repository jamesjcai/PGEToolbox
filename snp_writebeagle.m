function [status] = snp_writebeagle(geno, mark, filename)
%SNP_WRITEBEAGLE - saves as file of format BEAGLE 3.2
%snp_writebeagle(geno,mark,filename)

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
if nargin < 2
    mark = [];
end
%if (isempty(mark)), status=0; return; end

if (nargin < 3),
    [filename, pathname, filterindex] = uiputfile( ...
        {'*.beagle;*.bgl', 'BEAGLE Format Files (*.beagle, *.bgl)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Save as');
    if ~(filename), status = 0;
        return;
    end
    filename = [pathname, filename];
    if (filterindex == 1)
        if (isempty(find(filename == '.', 1)))
            filename = [filename, '.bgl'];
        end
        filenamemarker = [filename, '.markinfo'];
    end
else
    filenamemarker = [filename, '.markinfo'];
end

fid = fopen(filename, 'wt');
if (fid == -1),
    status = 0;
    warning('Unable to open file.');
    return;
end


[samplen, marklen] = snp_samplen(geno);
indvlen = samplen / 2;
transposedgeno = snp_genotranspose(geno);

ACGT = 'ACGT?';

for k = 1:marklen
    if ~isempty(mark)
        fprintf(fid, 'M %s ', mark.rsid{k});
    else
        fprintf(fid, 'M mark%d ', k);
    end

    for j = 1:size(transposedgeno, 2) - 1
        fprintf(fid, '%s ', ACGT(transposedgeno(k, j)));
    end
    fprintf(fid, '%s', ACGT(transposedgeno(k, end)));
    fprintf(fid, '\n');
end
fclose(fid);

%{
if (~isempty(mark))
    fid = fopen(filenamemarker,'wt');
    %fprintf(fid,'MARKER_ID\tSNP_rs#\tbp_POSITION\n');
    for k=1:marklen
        %fprintf(fid,'%d\t%s\t',k,mark.rsid{k});
        %fprintf(fid,'%d\t',k);
        fprintf(fid,'%s\t',mark.rsid{k});
        fprintf(fid,'%d\t',mark.pos(k));
        a=mark.allele{k};
        fprintf(fid,'%s\t%s\n',a(1),a(3));
    end
    fclose(fid);
end
%}

status = 1;
