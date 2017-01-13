function [s0,markinfo] = snp_downloadhapmapfreq(query,popcode,noise)
%SNP_DOWNLOADHAPMAPFREQ - downloads SNP frequency data from HapMap website
%[s0,genodata,markinfo] = snp_downloadhapmap(query,popcode,noise)

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

if (isempty(query)|isempty(popcode)), return; end

	popcode=upper(popcode);
	if ~(ismember(popcode,{'CEU','CHB','JPT','YRI','JPT+CHB'})),
	  error('Not available population code.');
	end
	if (strcmp(popcode,'JPT+CHB'))
	      popcode='JPT%2BCHB';
    end

base='http://www.hapmap.org/cgi-perl/gbrowse/hapmap24_B36/';
%base='http://www.hapmap.org/cgi-perl/gbrowse/gbrowse/hapmap/';
%base='http://www.hapmap.org/cgi-perl/gbrowse/hapmap20_B35/';
%base='http://www.hapmap.org/cgi-perl/gbrowse/hapmap_B34/';

base = strcat(base, '?plugin=SNPAlleleFrequencyDumper&plugin_action=Go&SNPAlleleFrequencyDumper.format=text&SNPAlleleFrequencyDumper.strand=fwd');
base = strcat(base, '&SNPAlleleFrequencyDumper.pop_code=', popcode);

%base = strcat(base, '?plugin=PhasedHaplotypeDumper&plugin_action=Go&PhasedHaplotypeDumper.format=text');
%base = strcat(base, '&PhasedHaplotypeDumper.pop_code=', popcode);

url = strcat(base, '&name=', query);



if (noise),
	disp('Downloading HapMap Genotype data from www.hapmap.org ... ')
	disp(url)
end

try
[s0,status]=urlread(url);
catch
	disp('Unable to download data.');
	return;
end

if ~(status>0),
	disp('Unable to download data.');
	return;
end


 	if (nargout>1)   % need genodata
		xx=strread(s0,'%s','delimiter','\n');
		for (k=4:length(xx)),
		      x=xx{k};

             try
		    [a,b]=regexp(x,'^rs\d+');
		    markinfo.rsid{k-3}=x(a:b);
		    [a,b]=regexp(x,'QC\+.+');
		    y=x(a+4:b);
		    %G 0.750 90 T 0.250 30 120
		    [a1,f1,n1,a2,f2,n2,t12]=strread(y,'%s%f%d%s%f%d%d');

		    if (f1>=f2)
			    markinfo.allele{k-3}=[a2{1},'/',a1{1}];
			    markinfo.maf(k-3)=f2;
             else
			    markinfo.allele{k-3}=[a1{1},'/',a2{1}];
			    markinfo.maf(k-3)=f1;
            end
             catch
                 markinfo=[];
                 return;
             end


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