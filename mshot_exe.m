function [a,pos]=mshot_exe(varstr,fileid)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if nargin<2
    fileid=1;
end
if nargin<1    
    %varstr='180 1 -s 388 -G 115.52 -eG 0.012 0.0 -eN 0.023 1'; 
    %varstr='180 1 -s 388 -G 115.52 -eG 0.012 0.0 -eN 0.023 1 -r 1000
    %750000';
    varstr='180 1 -s 388 -r 200 75000 -v 1 35000 40000 10';
end

pw0=pwd;
pw1=fileparts(which(mfilename));
cd(pw1);

filename=sprintf('output%d.txt',fileid);
if isunix
    cmdstr=sprintf('./msHOT %s >%s', varstr, filename);
elseif ispc
    cmdstr=sprintf('msHOT.exe %s >%s', varstr, filename);
elseif ismac
    cmdstr=sprintf('./msHOT %s >%s', varstr, filename);
else
    cmdstr=sprintf('./msHOT %s >%s', varstr, filename);
end
system(cmdstr);

[OUT]=readmsoutput(filename);
a=OUT.gametes{1};
pos=OUT.positions{1};
cd(pw0);