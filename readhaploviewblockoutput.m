function [data]=readhaploviewblockoutput(filename)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

    %txt = textread(filename,'%s','delimiter','\n','whitespace','');
    fid=fopen(filename);
    txt = textscan(fid,'%s','delimiter','\n','whitespace','','bufsize',4095*10);
    txt=txt{1};
    
    mt=cellfun('isempty',txt);
    txt(mt) = [];
    txt{end+1}='BLOCK';
    
    idx=strncmpi(txt,'Multiallelic',12);
    txt(idx)=[];
    
    idx=find(strncmpi(txt,'BLOCK',5));
    N=length(idx)-1;    % number of blocks
    %data{N}.markers = [];
    data=cell(1,N);    
    for k=1:length(idx)-1
        %h1=idx(k):idx(k+1);
        h1=txt{idx(k)};
        markers=sscanf(h1(strfind(h1,':')+1:end),'%d');
        data{k}.markers=uint32(markers');        
        nhap=idx(k+1)-idx(k)-1;
        nmak=size(markers,1);
        
        hapfreq=ones(1,nhap,'single');
        hapdata=ones(nhap,nmak,'uint8');
        ix=1;
        for kk=idx(k)+1:idx(k+1)-1
            %fprintf('%s\n',txt{kk})
            h2=txt{kk};
            pos1=strfind(h2,'(');            
            b=uint8(i_char2num(h2(1:pos1-2)));            
            hapfreq(1,ix)=single(str2double(h2(pos1+1:strfind(h2,')')-1)));
            hapdata(ix,:)=b;
            ix=ix+1;
        end
        data{k}.hapfreq=hapfreq;
        data{k}.hapdata=hapdata;
    end
    fclose(fid);
    
    
function [d]=i_char2num(s)
    [n]=length(s);
    d=zeros(1,n);
    for k=1:n
        d(k)=str2double(s(k));
    end
    