function [data]=snp_ldblock(genodata,gmarkinfo,method)
%SNP_LDBLOCK - returns haplotype block definitions (by Gabriel et al. method)
%
% Syntax: [data]=snp_ldblock(genodata,gmarkinfo);
%
%   S = SNP_LDBLOCK(G,M) reads a genotype data G and marker information M, 
%   returning the haplotype block data as a structure. S.markers is the marker 
%   index. S.hapdata is the haplotye data. S.hapfreq is the frequency of each
%   haplotype in the block.

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<3
    method='GAB';
end

    
% validmethods={'GAB','GAM','SPI'};
% for details, see the "Block Output Options" section of the Haploview
% Manual.

switch upper(method)
    case 'GAB'   % (Gabriel et al)
        outfile='input.ped.GABRIELblocks';
    case 'GAM'   % (4 gamete blocks)
        outfile='input.ped.4GAMblocks';
    case 'SPI'   % (solid spine blocks)
        outfile='input.ped.SPINEblocks';
    otherwise
        error('Invalid METHOD option.');
end
if ~isfield(gmarkinfo,'rsid')
    gmarkinfo.rsid=num2cellstr(1:length(gmarkinfo.pos));
end
        
%haploviewrun(genodata,gmarkinfo,method);

oldpath=pwd;
cdpge;
cd('addins/Haploview');

if exist('input.ped','file'), delete('input.ped'); end
snp_writelinkage(genodata,gmarkinfo,'input.ped');

%fid=fopen('input.map','w');
%for k=1:length(gmarkinfo.pos)
%    fprintf(fid,'%s\t%d\n',gmarkinfo.rsid{k},gmarkinfo.pos(k));
%end
%fclose(fid);

cmdline=sprintf('java -jar Haploview.jar -n -pedfile input.ped -info input.ped.map -skipcheck -blockoutput %s',method);
if exist(outfile,'file'), delete(outfile); end

[~,~]=system(cmdline);
[data]=i_parseblockfile(outfile,gmarkinfo);
cd(oldpath);


function [data]=i_parseblockfile(filename,gmarkinfo)
    %txt = textread(filename,'%s','delimiter','\n','whitespace','');
    fid = fopen(filename,'r');
    if fid==-1,
        data=[];
        return;
    end
    InputText=textscan(fid,'%s','delimiter','\n'); 
    txt=InputText{1};
    fclose(fid);

    mt=cellfun('isempty',txt);
    txt(mt) = [];
    txt{end+1}='BLOCK';
    
    idx=strncmpi(txt,'Multiallelic',12);
    txt(idx)=[];
    
    idx=find(strncmpi(txt,'BLOCK',5));
    N=length(idx)-1;    % number of blocks

    data(N,1).markers=[];
    
    for k=1:length(idx)-1
        %h1=idx(k):idx(k+1);
        h1=txt{idx(k)};
        markers=sscanf(h1(strfind(h1,':')+1:end),'%d');
        data(k,1).markers=uint16(markers');
        data(k,1).pos=uint32(gmarkinfo.pos(markers));
        hapfreq=[];
        hapdata=[];
        for kk=idx(k):idx(k+1)-1
            %fprintf('%s\n',txt{kk})
            h2=txt{kk};
            pos1=strfind(h2,'(');
            p=str2num(h2(pos1+1:strfind(h2,')')-1));
            b=i_char2num(h2(1:pos1-2));
            hapfreq=[hapfreq,p];
            hapdata=[hapdata;b];
        end
        data(k,1).hapfreq=single(hapfreq);
        data(k,1).hapdata=uint8(hapdata);
    end    
    
    
function [d]=i_char2num(s)
    [n]=length(s);
    d=zeros(1,n);
    for k=1:n
        d(k)=str2double(s(k));
    end
    