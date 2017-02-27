function haploviewrun(genodata,gmarkinfo,method)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-04-23 19:44:14 -0500 (Tue, 23 Apr 2013) $
% $LastChangedRevision: 529 $
% $LastChangedBy: jcai $

if nargin<3
    method='GAB';
end
    % validmethods={'GAB','GAM','SPI'};
    % for details, see the "Block Output Options" section of the Haploview
    % Manual.
    switch upper(method)
        case 'GAB'   % (Gabriel et al)
            outfile='input.ped.GABRIELblocks';
        case 'GAM'   % (4 gamete blocks)
            outfile='input.ped.4GAMblocks';
        case 'SPI'   % (solid spine blocks)
            outfile='input.ped.SPINEblocks';
        otherwise
            error('Invalid METHOD option.');
    end


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
cmdline=sprintf('java -jar Haploview.jar -n -pedfile input.ped -info input.map -skipcheck -blockoutput %s',method);
system(cmdline);

%[data]=i_parseblockfile(outfile);

cd(oldpath);

