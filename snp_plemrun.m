function [haplodata,markinfo]=snp_plemrun(genodata,fileidx,shownoise)
%SNP_PLEMRUN - runs PLEM
%
%[haplodata]=snp_plemrun(genodata)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-04-23 18:40:46 -0500 (Tue, 23 Apr 2013) $
% $LastChangedRevision: 526 $
% $LastChangedBy: jcai $

if nargin<3
    shownoise=true;
end
if nargin<2
    fileidx=1;
end
if nargin<1
    error('SNP_PLEMRUN:INPUT','Genotype input required.');
end

haplodata=[];
oldpath=pwd;
cdpge; cd('addins/plem');
%[exedir,dlgshown]=pge_getprgmdir(sprintf('%s_prgmdir',mfilename));
%if isempty(exedir)||dlgshown, return; end
%cd(exedir);

if shownoise, fprintf('Writing PLEM input file...'); end
inputfile=sprintf('input_%d.txt',fileidx);
outputfile=sprintf('output_%d.txt',fileidx);

[status]= snp_writeplem(genodata,inputfile);
if status==1    % good
    if shownoise, fprintf('done.\n'); end
else
    cd(oldpath);
    return;
end

if ispc
    cplotcmd=sprintf('plem.exe %s %s 0 2 50 20',inputfile,outputfile);
elseif ismac
    cplotcmd=sprintf('./plem.MAC %s %s 0 2 50 20',inputfile,outputfile);
else
    cplotcmd=sprintf('./plem.LINUX %s %s 0 2 50 20',inputfile,outputfile);
end

if shownoise, fprintf('Running PLEM ... '); end
try
    [~,~]=system(cplotcmd);
catch ME
    warning(ME.message);
    delete(inputfile); delete(outputfile);
    cd(oldpath);
    return;
end
if shownoise, fprintf('done.\n'); end
[~,majalle,minalle]=snp_maf(genodata);

if nargout>1
    [hap,markinfo]=snp_readplemout(outputfile);
    markinfo.majalle=uint8(majalle);
    markinfo.minalle=uint8(minalle);
elseif nargout==1
    [hap]=snp_readplemout(outputfile);
end

haplodata=uint8(zeros(size(hap)));
for k=1:size(hap,2)
    haplodata(hap(:,k),k)=minalle(k);
    haplodata(~hap(:,k),k)=majalle(k);
end
delete(inputfile);
delete(outputfile);
cd(oldpath);

