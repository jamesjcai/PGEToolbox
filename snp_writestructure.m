function [status]=snp_writestructure(geno,mark,filename)
%SNP_WRITESTRUCTURE - saves as STRUCTURE input file format
%
% snp_writestructure(geno,mark)
% snp_writestructure(geno,mark,filename)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (isempty(geno)), status=0; return; end
if (nargin < 3),
    [filename, pathname,filterindex] = uiputfile( ...
       {'*.structure;*.txt', 'STRUCTURE Input Files (*.structure, *.txt)';
        '*.*',  'All Files (*.*)'}, ...
        'Save as');
	if ~(filename), status=0; return; end
	filename=[pathname,filename];
	if (filterindex==1)
		if (isempty(find(filename=='.')))
			filename=[filename,'.txt'];
		end
	end
end
fid = fopen(filename,'wt');
if (fid == -1),
   status=0;
   warning('Unable to open file.');
   return;
end



[samplen,marklen]=snp_samplen(geno);
indvlen=samplen/2;

if (nargin < 2 || ~isfield(mark,'idvname')),
    for k=1:indvlen
      mark.idvname{k}=sprintf('Idv%d',k);
    end
end


%fprintf(fid,['%d\n'],indvlen);
%fprintf(fid,['%d\n'],marklen);
%fprintf(fid,['P %s\n'],sprintf('%d ',mark.pos));
%fprintf(fid,[char(ones(1,marklen)*['S']),'\n']);

%[geno] = snp_12geno(geno);
%{
%for (k=1:indvlen),
%      fprintf(fid,['#%d\n'],k);
%      for (j=1:marklen-1),
%        if (geno(k,j)==5)
%	      fprintf(fid,['?']);
%        else
%	      fprintf(fid,['%d'],geno(k,j));
%        end
%      end
%      fprintf(fid,['\n']);
%      for (j=1:marklen-1),
%        if (geno(k,j)==5)
%	    fprintf(fid,['?']);
%        else
%	    fprintf(fid,['%d'],geno(k,j)==0);
%        end
%      end
%	 fprintf(fid,['%d\n'],geno(k,j+1));
%end
%}


ACGT=[1 2 3 4 -9];
% -9 is undetermined
for (k=1:indvlen),
      fprintf(fid,'%s ', mark.idvname{k});
      for (j=1:2:marklen*2),
	      fprintf(fid,'%d ',ACGT(geno(k,j)));
      end
      fprintf(fid,'\n');

      fprintf(fid,'%s ', mark.idvname{k});
      for (j=2:2:marklen*2),
	      fprintf(fid,'%d ',ACGT(geno(k,j)));
      end
      fprintf(fid,'\n');
end

fclose(fid);
status=1;
