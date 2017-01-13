function snp_wirtetab(genodata,filename)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

fid = fopen(filename,'wt');
if (fid == -1), error('Unable to open file.'); end

[indvlen,marklen2]=size(genodata);
marklen=marklen2/2;

fprintf(fid,'%d\t%d\n',indvlen,marklen);

for (k=1:indvlen),
      fprintf(fid,'%d\t',k);
      for (j=1:marklen2),
	      fprintf(fid,'%d\t',genodata(k,j));
      end      
      fprintf(fid,'\n');
end
fclose(fid);

