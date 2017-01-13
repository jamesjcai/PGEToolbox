function snp_haploview(genodata,gmarkinfo)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-04-23 19:44:14 -0500 (Tue, 23 Apr 2013) $
% $LastChangedRevision: 529 $
% $LastChangedBy: jcai $

oldpath=pwd;
cdpge; cd('addins/Haploview');

%[exedir,dlgshown]=pge_getprgmdir(sprintf('%s_prgmdir',mfilename));
%if isempty(exedir)||dlgshown, return; end
%cd(exedir);

snp_writelinkage(genodata,gmarkinfo,'input.ped');
fid=fopen('input.map','w');
for k=1:length(gmarkinfo.pos)
    fprintf(fid,'%s\t%d\n',gmarkinfo.rsid{k},gmarkinfo.pos(k));
end
fclose(fid);
cmdline=sprintf('java -jar Haploview.jar -pedfile input.ped -info input.map -skipcheck');
system(cmdline);

%[data]=i_parseblockfile(outfile);

cd(oldpath);

