function cuffgtf2bed(gtf_file, fpkm_file, bed_file, color_value, trackname, description)
%CUFFGTF2BED Convert gtf to bed
%   cuffgtf2ucscbed(gtf_file, fpkm_file, bed_file, color_value) converts the gtf file generated by cufflink to bed file for UCSC genome browser
%           gtf_file            gtf file generated by the Cufflink.
%           fpkm_file           calculated FPKM value of isoforms.
%           bed_file            output file
%           color_value         optional, default setting is '0,0,0'. Color setting for the entire track in RGB form without internal space, e.g. '100,50,80', '255,255,255'.
%           trackname           optional, default setting is 'usertrackname'. Name displayed in UCSC genome browser custom track.
%           description           optional, default setting is 'user track description'. Description of custom track.
%   NOTE: other settings for bed file.
%          useScore=1         determine the displayed level of gray basing on score value
%          visibility=full    track displayed mode
%          score              score values are calculated according to FPKM value of each isoform. Internal boundary prefers the larger value, e.g. the score is 222 for FPKM 0.5.
%                             FPKM   0-0.5    0.5-1    1-2    2-4    4-8    8-16    16-32    32-64    64+
%                             score     111       222     333   444   555    666       777       888     999
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if isempty(color_value) color_value='0,0,0'; end
if isempty(trackname) trackname='usertrackname'; end
if isempty(description) description='user track description'; end

fid=fopen(fpkm_file);
content=textscan(fid,'%s%s%s%s%s%s%s%s%s%f%f%f%s','delimiter','\t','HeaderLines',1);
fclose (fid);
isoform=content{1};
fpkm=content{10};
gray_boundary=[0 0.5 1 2 4 8 16 32 64];
gray_value=[111 222 333 444 555 666 777 888 999];
score=zeros(1,length(fpkm));
for i=1:length(fpkm)    
    score(i)=gray_value(sum(fpkm(i)>=gray_boundary));
end



fid=fopen(gtf_file);
content=textscan(fid,'%s%s%s%s%s%s%s%s%s','delimiter','\t');
fclose (fid);
temp=regexp(content{9},'gene_id "([^"]+)"; transcript_id "([^"]+)";','tokens');
for i=1:length(temp)
    ID{i}=temp{i}{1}{2};
    Name{i}=temp{i}{1}{1};
end



fid=fopen(bed_file,'w');
fprintf(fid,'track name=%s description="%s" useScore=1 color=%s visibility=pack\n', trackname, description, color_value);

IDs=uunique(ID);
for i=1:length(IDs)
    tag_ID=IDs{i};
    % -------------------------------------
    list=find(strcmp(tag_ID,ID));
    %list=find(ismember(ID,tag_ID));
    % -------------------------------------    
    for j=1:8
        eval(sprintf('c%d=content{%d}(list);',j,j));
    end
    tag_name=unique(Name(list));
    tag_name=tag_name{1};
    
    if length(unique(c1))~=1
        c1
        error(sprintf('multiple chromosome ID in one transcription %s.', tag_ID));
    else
        c1=unique(c1);
        c1=c1{1};
    end

    if length(unique(c3))~=2
        unique(c3)
        
        error(sprintf('more type besides exon in %s.', tag_ID));
    end
    if length(unique(c7))~=1
        c7
        error(sprintf('multiple strand in one transcription %s.', tag_ID));
    else
        c7=unique(c7);
        c7=c7{1};
    end
    
    sublist=find(strcmp(c3,'exon'));
    temp4=c4(sublist);
    temp5=c5(sublist);
    exon=[];
    for j=1:length(temp4)
        exon=[exon; str2num(temp4{j}) str2num(temp5{j})];
    end
    if isempty(exon)
        error(sprintf('empty exon in transcription %s.', tag_ID));
    end
    exon=sort(exon,1);
    exon=sort(exon,2);
    blocksize=sprintf('%d',exon(1,2)-exon(1,1)+1);
    blockstart=sprintf('%d',exon(1,1)-exon(1));
    for j=2:size(exon,1)
        blocksize=[blocksize sprintf(',%d',exon(j,2)-exon(j,1)+1)];
        blockstart=[blockstart sprintf(',%d', exon(j,1)-exon(1))];
    end    
    fprintf(fid,'%s\t%d\t%d\t%s\t%d\t%s\t0\t0\t0\t%d\t%s\t%s\n',c1, min(min(exon))-1, max(max(exon)),...
        tag_ID, score(find(strcmp(tag_ID,isoform))), c7, size(exon,1), blocksize, blockstart);

end

fclose(fid);

end



function [unsorteduniques ia ib] = uunique(vec) 
    vec = vec(:)'; 
    [v a b] = unique(vec, 'first'); 
    if nargout > 2 
        [ia v] = sort(a); 
        [v ib] = ismember(b, v); 
    else 
       ia = sort(a); 
    end 
    unsorteduniques = vec(ia); 
end