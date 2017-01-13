function [haplodata,hmarkinfo]=snp_fastphaserun(genodata,gmarkinfo)
%SNP_FASTPHASERUN - runs FastPHASE
%
%[haplodata,hmarkinfo] = snp_fastphaserun(genodata,gmarkinfo)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2014-03-26 09:30:26 -0500 (Wed, 26 Mar 2014) $
% $LastChangedRevision: 758 $
% $LastChangedBy: jcai $

haplodata=[]; hmarkinfo=[];
oldpath=pwd;
cdpge(true); cd('addins/fastPhase');
%[exedir,dlgshown]=pge_getprgmdir(sprintf('%s_prgmdir',mfilename));
%if isempty(exedir)||dlgshown, return; end
%cd(exedir);


fprintf('Writing fastPHASE input file...');
[status]= snp_writephase(genodata,gmarkinfo,'input.inp');
if status==1    % good
    fprintf('done.\n');
else
    cd(oldpath);
    return;
end

if ispc
    cmdline='fastPHASE.exe input.inp output.txt';
elseif ismac
    cmdline='./fastphase_mac input.inp output.txt';
else
    cmdline='./fastphase_linux input.inp output.txt';
end


fprintf('Running fastPHASE ...');
try
    system(cmdline);
catch ME
    warning(ME.message,'fastPhase may not have been compiled properly on this machine.');
    cd(oldpath);
    return;
end
    fprintf('done.\n');
try
[haplodata,hmarkinfo] = snp_readfastphaseout('output.txt',1);
hmarkinfo.pos=gmarkinfo.pos;
hmarkinfo.rsid=gmarkinfo.rsid;
catch ME
    warning(ME.message,'Parsing addins/fastPhase/output.txt error.');
    cd(oldpath);
    return;
end
cd(oldpath);

