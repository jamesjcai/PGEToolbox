function [status] = snp_writemsout(genodata,ancalle,filename)
%snp_writemsout - save genotype data as MS output
% [status] = snp_writemsout(genodata,ancalle,filename)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

status=0;
[n,m]=size(genodata);
smpln=n*2;
markn=m/2;

%[smpln,markn]=snp_samplen(genodata);

if nargin < 3
    [filename, pathname,filterindex] = uiputfile( ...
       {'*.*'},'Save as');
	if ~(filename), status=0; return; end
	filename=[pathname,filename];
end

if nargin<2
    [p,majalle] = snp_maf(genodata);
    ancalle=majalle;
end

[genodata] = snp_01geno(genodata,ancalle);
[haplodata]=snp_geno2hap(genodata);


filename





fid1=fopen([filename,'.out'],'w');

%pos=markinfo.pos;
%rsid=markinfo.rsid;

%First line: any character. Use this line to store information about your
%data.

fprintf(fid1, '#NEXUS\n\n');
fprintf(fid1, '//\n');
fprintf(fid1, 'segsites:\n');
fprintf(fid1, 'positions:\n');
for k=1:size(haplodata,1)
    fprintf(fid1,'%d',haplodata(k,:));
    fprintf(fid1,'\n');
end

fclose(fid1);
status=1;