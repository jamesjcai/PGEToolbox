function [marker,popcode1,popcode2]=selectMarker2Popcodes(defaulttxt,idx1,idx2)

if nargin<3
    idx2=4;
end

if nargin<2
    idx1=1;
end

if nargin<1
    defaulttxt='';
end
items(1).name = 'MARKER or REGION:   ';
items(1).default = 0;        
items(1).values = {defaulttxt};
items(1).help = 'Search using a sequence name, gene name, locus, or other landmark.';


popset={'CEU','CHB','JPT','YRI','JPT+CHB'};
items(2).name = 'POPULATION 1:   ';
items(2).default = idx1;
items(2).indent = 0;
items(2).values = popset;
items(2).help = sprintf(['Population descriptors:\nCEU: CEPH (Utah residents with ancestry from northern and western Europe),\nCHB: Han Chinese in Beijing, China,\nJPT: Japanese in Tokyo, Japan,\nYRI: Yoruba in Ibadan, Nigeria.']);

items(3).name = 'POPULATION 2:   ';
items(3).default = idx2;
items(3).indent = 0;
items(3).values = popset;
items(3).help = sprintf(['Population descriptors:\nCEU: CEPH (Utah residents with ancestry from northern and western Europe),\nCHB: Han Chinese in Beijing, China,\nJPT: Japanese in Tokyo, Japan,\nYRI: Yoruba in Ibadan, Nigeria.']);


%title = 'Download from WWW.HAPMAP.ORG';
title = 'Download Data';
msg = sprintf(['The wildcard character * is allowed. \n\nExamples: Chr9:660000..760000, SNP:rs6870660,\nNM_153254, BRCA2, 5q31, ENm010.\n']);
%msg = 'xx';
out = CSEFlagDialog(items, title, msg);

if ~(isempty(out)),
	marker=deblank(out(1).answer);
	popcode1=popset{out(2).answer};
	popcode2=popset{out(3).answer};    
else
    marker='';
    popcode1='';
    popcode2='';    
end
