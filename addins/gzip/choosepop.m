function Index=choosepop(x)

filename='PopID.vcf';
fid=fopen(filename);
human_id={};
human_population={};
while ~feof(fid)
    temp=textscan(fgetl(fid),'%s','delimiter','\t');
    human_id=[human_id;temp{1}{1}];
    human_population=[human_population; temp{1}{2}];
end
races=unique(human_population);
fid=fopen('infile.vcf');
l=fgetl(fid);
while strncmp(l,'##',2)
    l=fgetl(fid);
end
tags=textscan(l,'%s');
tags=tags{1};
content=textscan(fid,repmat('%s',1,length(tags)),'delimiter','\t');
fclose(fid);
 
pop=x;
Index=[];
for i=1:length(pop)
    if ~ismember(pop{i},races) error('ERROR: unknown race'); end
        idx=strcmp(pop{i},human_population);
        [~,idx]=ismember(human_id(idx),tags);
        Index=[Index;idx];
end
Index=sort(Index)-9;
