function [genodata, markinfo, filename] = snp_downloadhapmap(query, popid, noise)
%SNP_DOWNLOADHAPMAP - downloads SNP genotype data from HapMap
%[s0,genodata,markinfo] = snp_downloadhapmap(query,popid,noise)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


genodata = [];
markinfo = [];
if (nargin < 3), noise = 0; end
if (nargin < 2),
    [query, popid] = selectMarkerPopcode;
end

if (isempty(query) || isempty(popid)), return; end

popid = upper(popid);
if ~(ismember(popid, {'CEU', 'CHB', 'JPT', 'YRI', 'JPT+CHB'})),
    error('Not available population code.');
end
if strcmp(popid, 'JPT+CHB')
    popid = 'JPT%2BCHB';
end

try
    base = urlread('http://www.bioinformatics.org/pgetoolbox/snp_downloadhapmap_url');
catch exception
    %warning('SNP_DOWNLOADHAPMAP can''t read online data. Default URL will be used.');
    disp(exception.message);
    base = 'http://hapmap.ncbi.nlm.nih.gov/cgi-perl/gbrowse/hapmap24_B36/';
end

base = strcat(base, '?plugin=SNPMultipopGenotypeDataDumper&plugin_action=Go&SNPMultipopGenotypeDataDumper.format=text');
base = strcat(base, '&SNPMultipopGenotypeDataDumper.pop_code=', popid);

url = strcat(base, '&name=', query);
if noise, fprintf('Downloading HapMap genotype data from URL:\n%s\n', url); end
fprintf('Downloading HapMap genotype data from URL:\n%s\n', url);

%{
try
    [s0,status]=urlread(url);
    if noise
        s0
    end
catch exception
    disp(exception.message);
    return;
end
if status~=1
    disp('Unable to download data.');
    return;
end
%}

if nargout % need genodata
    tempf = tempname;
    [filename, status] = urlwrite(url, tempf);
    pause(2);
    %[fid,Msg] = fopen(tempf,'wt');
    %if fid == -1, error(Msg); end
    %disp('Saving HapMap genotype data to a temporary file ... ')
    %fprintf(fid,'%s',s0);
    %fclose(fid);
    if status == 1
        [genodata, markinfo] = snp_readhapmap(filename);
    else
        disp('Unable to download data.');
        return;
    end
else
    [filename, pathname, filterindex] = uiputfile( ...
        {'*.hapmap;*.hmp', 'HapMap Genotype Files (*.hapmap, *.hmp)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Save as');
    if ~(filename), return; end
    filename = [pathname, filename];
    if filterindex == 1
        if isempty(find(filename == '.'))
            filename = [filename, '.hmp'];
        end
    end
    [~, status] = urlwrite(url, filename);
    %[fid,Msg] = fopen(filename,'wt');
    %if fid == -1, error(Msg); end
    %fprintf(fid,'%s',s0);
    %fclose(fid);
end

if strcmp(popid, 'JPT+CHB')
    markinfo.popid = 'JPT+CHB';
end