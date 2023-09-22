function [s0, genodata, markinfo] = snp_downloadperlegen(query, popcode, noise)
%SNP_DOWNLOADPERLEGEN - downloads SNP genotype data from Perlegen

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

s0 = '';
genodata = [];
markinfo = [];
if (nargin < 3), noise = 0; end
if (nargin < 2),
    [query, popcode] = selectMarkerPopcode;
end

if (isempty(query) | isempty(popcode)), return; end

popcode = upper(popcode);
if ~(ismember(popcode, {'CEU', 'CHB', 'JPT', 'YRI', 'JPT+CHB'})),
    error('Not available population code.');
end
if (strcmp(popcode, 'JPT+CHB'))
    popcode = 'JPT%2BCHB';
end

popcode = 'European+American';

%base='http://www.hapmap.org/cgi-perl/gbrowse/gbrowse/hapmap/';
%base='http://www.hapmap.org/cgi-perl/gbrowse/hapmap20_B35/';
base = 'http://genome.perlegen.com/cgi-bin/gbrowse/build34';
base = strcat(base, '?HaploviewDumper.disposition=view&plugin_action=Go&plugin_action=Go&plugin=HaploviewDumper');
base = strcat(base, '&HaploviewDumper.population=', popcode);
url = strcat(base, '&name=', query);


if (noise),
    disp('Downloading genotype data from genome.perlegen.com ... ')
    disp(url)
end


try
    [s0, status] = urlread(url);
catch
    disp('Unable to download data.');
    return;
end
if ~(status > 0),
    disp('Unable to download data.');
    return;
end

if (nargout > 1) % need genodata
    tempf = tempname;
    [fid, Msg] = fopen(tempf, 'wt');
    if fid == -1, error(Msg); end
    disp('Saving Perlegen data to a temporary file ... ')
    fprintf(fid, '%s', s0);
    fclose(fid);
    [genodata, markinfo] = snp_readperlegen(tempf);
end

if (nargout == 0) % don't need genodata, just save
    [filename, pathname, filterindex] = uiputfile( ...
        {'*.hapmap;*.hmp', 'HapMap Genotype Files (*.hapmap, *.hmp)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Save as');
    if ~(filename), return; end
    filename = [pathname, filename];
    if (filterindex == 1)
        if (isempty(find(filename == '.')))
            filename = [filename, '.hmp'];
        end
    end

    [fid, Msg] = fopen(filename, 'wt');
    if fid == -1, error(Msg); end
    fprintf(fid, '%s', s0);
    fclose(fid);
end