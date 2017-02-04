function [bp,LMap] = chrlen_hg19(id)
%CHRLEN - returns human chromosome length
%USAGE: [bp] = chrlen(id)
%
%-Version: Assembly Statistics for GRCh37 Release date: February 27, 2009

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if ischar(id)
	id=upper(id);
	if (strcmp(id,'23')), id='X'; end
	if (strcmp(id,'24')), id='Y'; end
    if (strcmp(id,'25')), id='MT'; end
elseif id==23
	id='X';
elseif id==24
	id='Y';
elseif id==24
	id='MT';
else
	id=int2str(id);
end

    
cnam = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14',...
        '15','16','17','18','19','20','21','22','X','Y','MT'};
    
clen=[249250621, 243199373, 198022430, 191154276, 180915260, 171115067,...
      159138663, 146364022, 141213431, 135534747, 135006516, 133851895,...
      115169878, 107349540, 102531392, 90354753, 81195210, 78077248,...
      59128983, 63025520, 48129895, 51304566, 155270560, 59373566, 16569];

[~,y]=ismember(upper(id),cnam);
if (y>0)
    bp=clen(y);
else
    error('%s is not a valid chromosome name.',id);
end


