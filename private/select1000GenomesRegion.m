function [region,chrid]=select1000GenomesRegion

chrset=textscan(num2str(1:22),'%s');
chrset=chrset{1};
chrset{23}='X';

items(1).name = 'CHROMOSOME:       ';
items(1).default = 1;
items(1).indent = 0;
items(1).values = chrset;
items(1).help = sprintf('Population descriptors:\nCEU: CEPH (Utah residents with ancestry from northern and western Europe),\nCHB: Han Chinese in Beijing, China,\nJPT: Japanese in Tokyo, Japan,\nYRI: Yoruba in Ibadan, Nigeria.');


items(2).name = 'START POSITION:   ';
items(2).default = 0;        
items(2).values = {'1'};
items(2).help = 'Search using a sequence name, gene name, locus, or other landmark.';

items(3).name = 'END POSITION:     ';
items(3).default = 0;        
items(3).values = {'50000'};
items(3).help = '(max 10K region)';



%title = 'Download from WWW.HAPMAP.ORG';
title = 'Data Slicer';
msg = sprintf('The Data Slicer provides an interface which allows users to get subsections of either vcf (vcftools) or bam (samtools) files based on genomic coordinates.\n\n\n');
out = CSEFlagDialog(items, title, msg);

if ~(isempty(out)),
	chrid=chrset{out(1).answer};
    startpos=str2double(out(2).answer);
    endpos=str2double(out(3).answer);
    %if endpos-startpos+1>10000
    %    errordlg('Max 10K region');
    %    chrid='';
    %    region='';
    %    return;
    %end
    region=sprintf('%s:%s-%s',chrid,...
        deblank(out(2).answer),deblank(out(3).answer));	
else
    region='';
    chrid='';
end
