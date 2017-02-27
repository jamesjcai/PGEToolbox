function [region,popcode]=selectHapMapRegion(defaulttxt)

warning('Example inputs are used.');
region='BRCA2'
popcode='CEU'


%{
if nargin<1
    defaulttxt='';
end
items(1).name = 'MARKER or REGION:   ';
items(1).default = 0;        
items(1).values = {defaulttxt};
items(1).help = 'Search using a sequence name, gene name, locus, or other landmark.';


popset={'CEU','CHB','JPT','YRI','JPT+CHB'};
items(2).name = 'POPULATION:   ';
items(2).default = 1;
items(2).indent = 0;
items(2).values = popset;
items(2).help = sprintf(['Population descriptors:\nCEU: CEPH (Utah residents with ancestry from northern and western Europe),\nCHB: Han Chinese in Beijing, China,\nJPT: Japanese in Tokyo, Japan,\nYRI: Yoruba in Ibadan, Nigeria.']);


%title = 'Download from WWW.HAPMAP.ORG';
title = 'Download Data';
msg = sprintf(['The wildcard character * is allowed. \n\nExamples: Chr9:660000..760000, SNP:rs6870660,\nNM_153254, BRCA2, 5q31, ENm010.\n']);
out = CSEFlagDialog(items, title, msg);

if ~(isempty(out)),
	region=deblank(out(1).answer);
	popcode=popset{out(2).answer};
else
    region='';
    popcode='';
end
%}