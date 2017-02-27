function [bp,LMap] = chrlen_hg18(id)
%CHRLEN - returns human chromosome length
%USAGE: [bp] = chrlen(id)
%
%-Version: from UCSC Genome Browser based on Human Mar. 2006 Assembly.
% http://genome.ucsc.edu/cgi-bin/hgTracks?hgsid=91767390&chromInfoPage=
% Human Mar. 2006 (hg18) Browser Sequences

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
elseif id==23
	id='X';
elseif id==24
	id='Y';
else
	id=int2str(id);
end

cnam = {'1','10','11','12','13','14','15','16','17','18','19','2','20',...
        '21','22','3','4','5','6','7','8','9','X','Y'};
clen=[247249719, 135374737, 134452384, 132349534, 114142980, 106368585,...
      100338915, 88827254, 78774742, 76117153, 63811651, 242951149,...
      62435964, 46944323, 49691432, 199501827, 191273063, 180857866,...
      170899992, 158821424, 146274826, 140273252, 154913754, 57772954];


[x,y]=ismember(upper(id),cnam);
if x
    bp=clen(y);
else
    error('%s is not a valid chromosome name.',id);
end


if nargout>1
	LMap={
'1',247249719;
'10',135374737;
'11',134452384;
'12',132349534;
'13',114142980;
'14',106368585;
'15',100338915;
'16',88827254;
'17',78774742;
'18',76117153;
'19',63811651;
'2',242951149;
'20',62435964;
'21',46944323;
'22',49691432;
'3',199501827;
'4',191273063;
'5',180857866;
'6',170899992;
'7',158821424;
'8',146274826;
'9',140273252;
'X',154913754;
'Y',57772954};
end





% #chrom	size	fileName
% chr1	247249719	/gbdb/hg18/nib/chr1.nib
% chr1_random	1663265	/gbdb/hg18/nib/chr1_random.nib
% chr10	135374737	/gbdb/hg18/nib/chr10.nib
% chr10_random	113275	/gbdb/hg18/nib/chr10_random.nib
% chr11	134452384	/gbdb/hg18/nib/chr11.nib
% chr11_random	215294	/gbdb/hg18/nib/chr11_random.nib
% chr12	132349534	/gbdb/hg18/nib/chr12.nib
% chr13	114142980	/gbdb/hg18/nib/chr13.nib
% chr13_random	186858	/gbdb/hg18/nib/chr13_random.nib
% chr14	106368585	/gbdb/hg18/nib/chr14.nib
% chr15	100338915	/gbdb/hg18/nib/chr15.nib
% chr15_random	784346	/gbdb/hg18/nib/chr15_random.nib
% chr16	88827254	/gbdb/hg18/nib/chr16.nib
% chr16_random	105485	/gbdb/hg18/nib/chr16_random.nib
% chr17	78774742	/gbdb/hg18/nib/chr17.nib
% chr17_random	2617613	/gbdb/hg18/nib/chr17_random.nib
% chr18	76117153	/gbdb/hg18/nib/chr18.nib
% chr18_random	4262	/gbdb/hg18/nib/chr18_random.nib
% chr19	63811651	/gbdb/hg18/nib/chr19.nib
% chr19_random	301858	/gbdb/hg18/nib/chr19_random.nib
% chr2	242951149	/gbdb/hg18/nib/chr2.nib
% chr2_random	185571	/gbdb/hg18/nib/chr2_random.nib
% chr20	62435964	/gbdb/hg18/nib/chr20.nib
% chr21	46944323	/gbdb/hg18/nib/chr21.nib
% chr21_random	1679693	/gbdb/hg18/nib/chr21_random.nib
% chr22	49691432	/gbdb/hg18/nib/chr22.nib
% chr22_random	257318	/gbdb/hg18/nib/chr22_random.nib
% chr22_h2_hap1	63661	/gbdb/hg18/nib/chr22_h2_hap1.nib
% chr3	199501827	/gbdb/hg18/nib/chr3.nib
% chr3_random	749256	/gbdb/hg18/nib/chr3_random.nib
% chr4	191273063	/gbdb/hg18/nib/chr4.nib
% chr4_random	842648	/gbdb/hg18/nib/chr4_random.nib
% chr5	180857866	/gbdb/hg18/nib/chr5.nib
% chr5_random	143687	/gbdb/hg18/nib/chr5_random.nib
% chr5_h2_hap1	1794870	/gbdb/hg18/nib/chr5_h2_hap1.nib
% chr6	170899992	/gbdb/hg18/nib/chr6.nib
% chr6_random	1875562	/gbdb/hg18/nib/chr6_random.nib
% chr6_cox_hap1	4731698	/gbdb/hg18/nib/chr6_cox_hap1.nib
% chr6_qbl_hap2	4565931	/gbdb/hg18/nib/chr6_qbl_hap2.nib
% chr7	158821424	/gbdb/hg18/nib/chr7.nib
% chr7_random	549659	/gbdb/hg18/nib/chr7_random.nib
% chr8	146274826	/gbdb/hg18/nib/chr8.nib
% chr8_random	943810	/gbdb/hg18/nib/chr8_random.nib
% chr9	140273252	/gbdb/hg18/nib/chr9.nib
% chr9_random	1146434	/gbdb/hg18/nib/chr9_random.nib
% chrM	16571	/gbdb/hg18/nib/chrM.nib
% chrX	154913754	/gbdb/hg18/nib/chrX.nib
% chrX_random	1719168	/gbdb/hg18/nib/chrX_random.nib
% chrY	57772954	/gbdb/hg18/nib/chrY.nib
