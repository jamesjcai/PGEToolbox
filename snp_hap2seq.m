function [seq] = snp_hap2seq(haplodata,hmarkinfo)
%SNP_HAP2SEQ - converts HAPLODATA to SEQ

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


pos=hmarkinfo.pos;
pos=pos-min(pos)+1;
[n]=size(haplodata,1);
seq=ones(n,max(pos)-min(pos)+1);
seq(:,pos)=haplodata;
