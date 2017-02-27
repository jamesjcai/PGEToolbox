function [status] = snp_writesweep(hapldata,markinfo,filename)
%HAP_WRITESWEEP - save haplotype data for a software called SWEEP
%
% hap_writesweep(hapldata,markinfo)
% hap_writesweep(hapldata,markinfo,filename)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

status=0;
if nargin < 3
    [filename, pathname,filterindex] = uiputfile( ...
       {'*.*'},'Save as');
	if ~(filename), status=0; return; end
	filename=[pathname,filename];
end

fid1=fopen([filename,'.snp'],'w');
fid2=fopen([filename,'.phase'],'w');

pos=markinfo.pos;
rsid=markinfo.rsid;
fprintf(fid1, 'snpid	chr	HG16\n');
warning('SNP_SAVE4SWEEP assumes all markers from chromosome 1.')
for k=1:length(rsid)
    fprintf(fid1, '%s\t1\t%d\n',rsid{k},pos(k));
end

[n,m]=size(hapldata);
for k=1:n
    fprintf(fid2, 'idv_%d\tT\t',k);
    for j=1:m-1
         fprintf(fid2,'%d\t',hapldata(k,j));
    end
    fprintf(fid2, '%d\n',hapldata(k,m));
end

fclose(fid1);
fclose(fid2);
status=1;
