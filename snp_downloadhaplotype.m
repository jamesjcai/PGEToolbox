function [hapldata,markinfo,filename] = snp_downloadhaplotype(query,popcode,noise)
%SNP_DOWNLOADHAPLOTYPE - downloads haplotype data from HapMap
%[s0,hapldata,markinfo] = snp_downloadhaplotype(query,popcode)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


hapldata=[]; markinfo=[];
if (nargin<3), noise=false; end
if (nargin<2),
	[query,popcode]=selectMarkerPopcode('Chr9:660000..760000');
end

if (isempty(query)||isempty(popcode)), return; end

	popcode=upper(popcode);
    if ~(ismember(popcode,{'CEU','CHB','JPT','YRI','JPT+CHB','JPT%2BCHB'})),
	  error('Not available population code.');
	end
	if (strcmp(popcode,'JPT+CHB'))
	      popcode='JPT%2BCHB';
    end

try
	base = urlread('http://www.bioinformatics.org/pgetoolbox/snp_downloadhapmap_url');
    %base ='http://hapmap.ncbi.nlm.nih.gov/cgi-perl/gbrowse/hapmap24_B36/';
    %base ='http://hapmap.ncbi.nlm.nih.gov/cgi-perl/gbrowse/hapmap3r3_B36/';
catch exception
    % warning('SNP_DOWNLOADHAPLOTYPE can''t read online data. Default URL will be used.');
    disp(exception.message);
    base='http://hapmap.ncbi.nlm.nih.gov/cgi-perl/gbrowse/hapmap24_B36/';
end

lowpopcode=lower(popcode);

if strcmp(popcode,'JPT')||strcmp(popcode,'CHB')
    popcode='JPT%2BCHB';
    lowpopcode='jc';
    disp('using JPT+CHB')
end

base = strcat(base, '?plugin=PhasedHaplotypeDumper&plugin_action=Go&PhasedHaplotypeDumper.format=text');
base = strcat(base, '&PhasedHaplotypeDumper.pop_code=', popcode);

%cons -+ Consensus Set (SNP in all three panels)
%poly -+ Polymorphic in at least one panel
%ceu - Polymorphic in CEU
%yri - Polymorphic in YRI
%jc - Polymorphic in JPT+CHB

%base = strcat(base, '&PhasedHaplotypeDumper.filters=', 'cons');
base = strcat(base, '&PhasedHaplotypeDumper.filters=', 'poly');
base = strcat(base, '&PhasedHaplotypeDumper.filters=', lowpopcode);

url = strcat(base, '&name=', query);
if noise, fprintf('Downloading HapMap haplotype data from URL:\n%s\n',url); end

%{
try
pause(1);    
[s0,status]=urlread(url);
catch
	disp('Unable to download data.');
	return;
end
if ~(status>0),
	disp('Unable to download data.');
	return;
end
%}
 	if nargout   % need hapldata
 	       tempf = tempname;
           pause(1)
           [filename,status]=urlwrite(url,tempf);
           if status==1
                [hapldata,markinfo] = snp_readhaplotype(filename);
           else
                disp('Unable to download data.');
                return;               
           end
    else
        [filename, pathname,filterindex] = uiputfile( ...
	       {'*.phased;*.phs', 'Phased Haplotype Files (*.phased, *.phs)';
		'*.*',  'All Files (*.*)'}, ...
		'Save as');
		if ~(filename), return; end
		filename=[pathname,filename];
		if (filterindex==1)
			if (isempty(find(filename=='.')))
			filename=[filename,'.phs'];
			end
        end
        [~,status]=urlwrite(url,filename);
    end
