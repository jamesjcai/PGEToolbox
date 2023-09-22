function [chrid, startn, endn] = chrloc2num(loctxt)
%chrloc2num - converts 'chr11:79399570-79399751' to 11,79399570,79399751

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-04-26 23:08:19 -0500 (Fri, 26 Apr 2013) $
% $LastChangedRevision: 532 $
% $LastChangedBy: jcai $

%loctxt='chr11:79399570-79399751';
loctxt = strrep(loctxt, ',', '');
loctxt = lower(loctxt);
idx = strfind(loctxt, ':');
chrid = str2double(loctxt(4:idx-1));

loctxt = loctxt(idx+1:end);
idx = strfind(loctxt, '-');
startn = str2double(loctxt(1:idx-1));
endn = str2double(loctxt(idx+1:end));
