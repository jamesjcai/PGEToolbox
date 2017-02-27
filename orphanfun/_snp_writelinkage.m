function [status] =  snp_writelinkage(geno,mark,filename)
%SNP_WRITELINKAGE - saves as linkage format
%snp_writelinkage(geno,mark,filename)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2012-12-09 22:08:33 -0600 (Sun, 09 Dec 2012) $
% $LastChangedRevision: 276 $
% $LastChangedBy: jcai $

if (isempty(geno)), status=0; return; end
if nargin<2
    mark=[];
end
%if (isempty(mark)), status=0; return; end

if (nargin < 3),
    [filename, pathname,filterindex] = uiputfile( ...
        {'*.pedigree;*.ped', 'Linkage Format Files (*.pedigree, *.ped)';
        '*.*',  'All Files (*.*)'}, ...
        'Save as');
	if ~(filename), status=0; return; end
	filename=[pathname,filename];
	if (filterindex==1)
		if (isempty(find(filename=='.')))
			filename=[filename,'.ped'];
        end
		filenamemarker=[filename,'.markinfo'];
    end
else
    filenamemarker=[filename,'.markinfo'];
end

fid = fopen(filename,'wt');
if (fid == -1),
   status=0;
   warning('Unable to open file.');
   return;
end

[samplen,marklen]=snp_samplen(geno);
indvlen=samplen/2;

ACGT='12340';

for (k=1:indvlen),
    
%fprintf(fid,['%d\n'],indvlen);
%fprintf(fid,['%d\n'],marklen);
%fprintf(fid,['P %s\n'],sprintf('%d ',mark.pos));
%fprintf(fid,[char(ones(1,marklen)*['S']),'\n']);
    
      fprintf(fid,['Family%d\tIndividual%d\t0\t0\t1\t0\t'],...
          k,k);
      
      for (j=1:marklen*2),
	      fprintf(fid,'%s\t',ACGT(geno(k,j)));
      end      
      fprintf(fid,'\n');
end
fclose(fid);

if (~isempty(mark))
fid = fopen(filenamemarker,'wt');
    %fprintf(fid,'MARKER_ID\tSNP_rs#\tbp_POSITION\n');
for k=1:marklen
    %fprintf(fid,'%d\t%s\t',k,mark.rsid{k});
    %fprintf(fid,'%d\t',k);
    fprintf(fid,'%s\t',mark.rsid{k});
    fprintf(fid,'%d\n',mark.pos(k));
end
fclose(fid);
end
status=1;


