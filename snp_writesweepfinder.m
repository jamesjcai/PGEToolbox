function [status] = snp_writesweepfinder(geno, mark, filename, ancalle)
%SNP_WRITESWEEPFINDER - saves GENODATA into SWEEPFINDER input file
%
% snp_writesweepfinder(geno,mark)
% snp_writesweepfinder(geno,mark,filename)
%
% SWEEPFINDER input file format
% ----------------------------------------------------------------
% The snp file should be a tab-delimited file with column headers, and
% one row per SNP.  One column header should be "x" (the frequency of the
% SNP), another should be "n" (the sample size, must be greater than x), and
% another should be "position" (the chromosomal location of the SNP).
% Optionally, a fourth column named "folded" can be added.  If it is present,
% than a value of one indicates that the SNP is folded (there is no distinction
% between ancestral/derived states), and 0 means unfolded.  If the folded
% column is not present, all SNPs are assumed to be unfolded.
% Column names do not actually contain quotes.  A sample input file might look
% something like:
%
% position        x       n	folded
% 37.000000       10      46	0
% 145.000000      3       47	0
% 277.000000      1       47	1
% 385.000000      37      43	1
% 469.000000      2       45	0
% 585.000000      1       44	0
% 733.000000      10      45	0
%
%SweepFinder program (Last updated: 8/14/06) Melissa Hubisz and Rasmus Nielsen
%Refererence: Nielsen et al, Genome Research 2005 Nov; 15(11):1566-75

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin < 4
    ancalle = [];
end
if (isempty(geno)), status = 0;
    return;
end
if (isempty(mark)), status = 0;
    return;
end
if (nargin < 3),
    [filename, pathname, filterindex] = uiputfile( ...
        {'*.tab;*.txt', 'SWEEPFINDER Input Files (*.tab, *.txt)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Save as');
    if ~(filename), status = 0;
        return;
    end
    filename = [pathname, filename];
    if (filterindex == 1)
        if (isempty(find(filename == '.', 1)))
            filename = [filename, '.tab'];
        end
    end
end
fid = fopen(filename, 'wt');
if (fid == -1),
    status = 0;
    warning('Unable to open file.');
    return;
end

[samplen, marklen] = snp_samplen(geno);
fprintf(fid, 'position\tx\tn\tfolded\n');
if isempty(ancalle)
    [maf, majalle, minalle] = snp_maf(geno);
    for k = 1:marklen
        if maf(k) > 0
            x = geno(:, [k * 2 - 1, k * 2]);
            x = x(:);
            n1 = sum(x == majalle(k));
            n2 = sum(x == minalle(k));
            fprintf(fid, '%d\t%d\t%d\t1\n', mark.pos(k), n2, n1+n2);
        end
    end
else
    error('Method with this option is under development.')
end
fclose(fid);
status = 1;
