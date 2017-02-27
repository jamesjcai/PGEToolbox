function [ancalle,chrmm]=snp_ancallechimp(chrid,pos,chrmm)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

ancalle=[];
if chrid>22 || chrid<1
    return;
end

if nargin<3
    chrmm=[];
chrmm = memmapfile(sprintf('Y:/Alignments/chimpaln/mmFilenamechr%d',chrid),...
    'format', 'uint8');
end
ancalle=chrmm.Data(pos);
%if ancalle==15
%    ancalle=5;
%end