function [rsid, mat] = snp_dbsnpquery(qstr)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

rsid = [];
mat = [];

if (nargin < 1),
    prompt = {sprintf('This command will return a list of SNPs in dbSNP.\n Please enter a search query: ')};
    def = {'((("8"[CHR]+AND+(5480000[Base+Position]+:+5500000[Base+Position]))+AND+"homo+sapiens"[Organism])+AND+"snp"[Snp_Class])'};
    dlgTitle = 'Search';
    lineNo = 1;
    answer = inputdlg(prompt, dlgTitle, lineNo, def);

    if ~(isempty(answer)),
        qstr = answer{1};
    else
        return;
    end
end

urlFetch = sprintf('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=snp&retmax=10&retmode=txt&term=%s', ...
    qstr);
try
    fetchResults = char(strread(urlread(urlFetch), '%s', 'delimiter', '\n', 'whitespace', ''));
    fetchResults = cellstr(fetchResults);
    fetchResults = strtrim(fetchResults);
    %fetchResults = deblank(fetchResults);
catch
    %errordlg(lasterr)
    disp(urlFetch)
    error(lasterr);
end


numLines = strmatch('<Count>', fetchResults);
[mat, idx] = regexp(fetchResults(numLines(1), :), '\d', 'match', 'start');
for i = 1:length(mat)
    mat{i} = char(mat{i})';
end
if str2num(mat{1}) > 1000
    warning('Over the maximum number (1000) of SNPs retriveled')
end


%select the ids from the xml
numLines = strmatch('<Id>', fetchResults);
[mat, idx] = regexp(fetchResults(numLines, :), '\d', 'match', 'start');

for i = 1:length(mat)
    mat{i} = char(mat{i})';
    rsid = [rsid, str2num(mat{i})];
end
