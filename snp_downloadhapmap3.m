function [s0,genodata,markinfo] = snp_downloadhapmap3(query,popcode,noise)
%SNP_DOWNLOADHAPMAP3 - downloads genotype data from HapMap3
%[s0,genodata,markinfo] = snp_downloadhapmap3(query,popcode,noise)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


s0=''; genodata=[]; markinfo=[];
if (nargin<3), noise=0; end
if (nargin<2),
	[query,popcode]=selectMarkerPopcode;
end

if (isempty(query)||isempty(popcode)), return; end

	popcode=upper(popcode);
	if ~(ismember(popcode,{'ASW','CEU','CHB','CHD','GIH','JPT','LWK','MEX','MKK','TSI','YRI'})),
	  error('Not available population code.');
    end
    
try
	%base = urlread('http://www.bioinformatics.org/pgetoolbox/snp_downloadhapmap3_url');
    base = 'http://hapmap.ncbi.nlm.nih.gov/cgi-perl/gbrowse/hapmap28_B36/';
catch
	warning('SNP_DOWNLOADHAPMAP3 can''t read online data. Default URL will be used.');
    base='http://hapmap.ncbi.nlm.nih.gov/cgi-perl/gbrowse/hapmap28_B36/';
end

%base='http://hapmap.ncbi.nlm.nih.gov/cgi-perl/gbrowse/hapmap27_B36/';
%base='http://www.hapmap.org/cgi-perl/gbrowse/hapmap27_B36/';
%base='http://www.hapmap.org/cgi-perl/gbrowse/hapmap3_B36/';
%base='http://www.hapmap.org/cgi-perl/gbrowse/gbrowse/hapmap/';
%base='http://www.hapmap.org/cgi-perl/gbrowse/hapmap20_B35/';
%base='http://www.hapmap.org/cgi-perl/gbrowse/hapmap_B34/';

base = strcat(base, '?plugin=SNPGenotypeDataPhase3Dumper&plugin_action=Go&SNPGenotypeDataPhase3Dumper.format=text');
base = strcat(base, '&SNPGenotypeDataPhase3Dumper.pop_code=', popcode);

%base = strcat(base, '?plugin=PhasedHaplotypeDumper&plugin_action=Go&PhasedHaplotypeDumper.format=text');
%base = strcat(base, '&PhasedHaplotypeDumper.pop_code=', popcode);

url = strcat(base, '&name=', query);

if (noise),
	disp('Downloading HapMap3 Genotype data from www.hapmap.org ... ')
	disp(url)
end

try
[s0,status]=urlread(url);

 tempf = tempname;
[filename,status]=urlwrite(url,tempf);
catch
	disp('Unable to download data.');
	return;
end

if ~(status>0),
	disp('Unable to download data.');
	return;
end

 	if (nargout>1)   % need genodata
 	       tempf = tempname;           
 	       [fid,Msg] = fopen(tempf,'wt');
 	       if fid == -1, error(Msg); end
           if noise
    	       disp('Saving HapMap genotype data to a temporary file ... ')
           end
 	       fprintf(fid,'%s',s0);
 	       fclose(fid);
           try
  	       [genodata,markinfo] = snp_readhapmap(tempf);
           catch
           [genodata,markinfo] = snp_readhapmap(filename);
           end
           
           if (~isempty(strfind(markinfo.popid,'Unknown')))
               markinfo.popid=popcode;
           end
	end

	if (nargout==0)   % don't need genodata, just save
	    [filename, pathname,filterindex] = uiputfile( ...
	       {'*.hapmap;*.hmp', 'HapMap Genotype Files (*.hapmap, *.hmp)';
		'*.*',  'All Files (*.*)'}, ...
		'Save as');
		if ~(filename), return; end
		filename=[pathname,filename];
		if (filterindex==1)
			if (isempty(strfind(filename,'.')))
			filename=[filename,'.hmp'];
			end
		end

	       [fid,Msg] = fopen(filename,'wt');
	       if fid == -1, error(Msg); end
	       fprintf(fid,'%s',s0);
	       fclose(fid);
    end