function [haplodata,hmarkinfo]=snp_phaserun(genodata,gmarkinfo)
%PHASERUN - runs PHASE
%
%[haplodata,hmarkinfo] = snp_phaserun(genodata,gmarkinfo)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2014-03-26 09:30:26 -0500 (Wed, 26 Mar 2014) $
% $LastChangedRevision: 758 $
% $LastChangedBy: jcai $

haplodata=[]; hmarkinfo=[];
oldpath=pwd;
cdpge(true); cd('addins/phase');
%[exedir,dlgshown]=pge_getprgmdir(sprintf('%s_prgmdir',mfilename));
%if isempty(exedir)||dlgshown, return; end
%cd(exedir);

if nargin<2
    gmarkinfo=[];
end

%fprintf('Writing PHASE input file...');
[status]= snp_writephase(genodata,gmarkinfo,'input.inp');
if status==1    % good
    %fprintf('done.\n');
else
    cd(oldpath);
    return;
end

if ispc
    cplotcmd='phase.exe input.inp output.txt';
elseif ismac
    cplotcmd='./phase_mac input.inp output.txt';
else
    cplotcmd='./phase_linux input.inp output.txt';
end


%fprintf('Running PHASE ... ');
try
[~,~]=system(cplotcmd);
catch ME
    warning('PHASE may not have been compiled properly on this machine.');
    cd(oldpath);
    return;
end
    %fprintf('done.\n');

[haplodata,hmarkinfo] = snp_readphaseout('output.txt',1);

hmarkinfo.pos=gmarkinfo.pos;
hmarkinfo.rsid=gmarkinfo.rsid;
%hmarkinfo.maf=gmarkinfo.maf;
cd(oldpath);

