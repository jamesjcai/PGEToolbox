function [a,pos]=ms_exe(varstr,fileid)

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
    %varstr='180 1 -s 388 -G 115.52 -eG 0.012 0.0 -eN 0.023 1 -r 1000 750000';     
    
    %varstr='68 1 -t 1882 -r 3764 50000 -I 2 24 44 160 -g 1 4303.90 -g 2 41228.0 -eg 0.001070 1 0. -eg 0.0001117 2 0. -en 0.000775 2 9.1e-05 -en 0.000785 2 0.01 -ej 0.001600 2 1';
    %varstr='68 1 -t 1792 -r 3584 50000 -I 2 24 44 160 -g 1 3704.88 -g 2 71464.5 -eg 0.001243 1 0. -eg 6.444e-05 2 0. -es 0.000590 2 0.95    -en 0.000590 3 0.01 -en 0.000775 2 7.14286e-05 -en 0.000785 2 0.01 -ej 0.0016 2 1 -ej 0.005 3 1';
    %REF: Plagnol V, Wall JD. PLoS Genet. 2006 Jul;2(7):e105.   
    varstr='simple';
end

pw0=pwd;
pw1=fileparts(which(mfilename));
cd(pw1);
cd('addins/ms_exe');

switch lower(varstr)
    case 'simple'
        varstr='120 1 -t 20';
   case {'exponential','exp'}        
        varstr='120 1 -t 6.4 -G 6.93 -eG 0.2 0.0 -eN 0.3 0.5';
    case 'bottleneck'
        varstr='120 1 -t 6.4 -G 6.93 -eG 0.2 0.0 -eN 0.3 0.5';
    case 'structure'
        varstr='120 1 -t 1882 -r 3764 50000 -I 2 60 60 1 -g 1 4303.90 -g 2 41228.0 -eg 0.001070 1 0. -eg 0.0001117 2 0. -en 0.000775 2 9.1e-05 -en 0.000785 2 0.01 -ej 0.001600 2 1';
    case 'wall'
        varstr='68 1 -t 1882 -r 3764 50000 -I 2 24 44 160 -g 1 4303.90 -g 2 41228.0 -eg 0.001070 1 0. -eg 0.0001117 2 0. -en 0.000775 2 9.1e-05 -en 0.000785 2 0.01 -ej 0.001600 2 1';
end

filename=sprintf('output%d.txt',fileid);
if isunix
    cmdstr=sprintf('./ms %s >%s', varstr, filename);
elseif ispc
    cmdstr=sprintf('ms.exe %s >%s', varstr, filename);
elseif ismac
    cmdstr=sprintf('./ms %s >%s', varstr, filename);
else
    cmdstr=sprintf('./ms %s >%s', varstr, filename);
end
system(cmdstr);


[OUT]=readmsoutput(filename);
a=OUT.gametes{1};
pos=OUT.positions{1};
delete(filename);
cd(pw0);


