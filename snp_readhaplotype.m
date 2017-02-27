function [hapldata,markinfo] = snp_readhaplotype(filename,noise)
%snp_readhaplotype - Reads HapMap haplotype data
%[genodata,markinfo] = snp_readhaplotype(filename,noise)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin < 1
    [filename, pathname] = uigetfile( ...
       {'*.phased;*.phs', 'HapMap Files (*.phased, *.phs)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a HapMap (phased haplotype data) format file');
	if ~(filename), hapldata=[]; markinfo=[]; return; end
	filename=[pathname,filename];
end

if nargin < 2, noise=0; end

if ~(exist(filename,'file') || exist(fullfile(cd,filename),'file')),
    error('InvalidInput','Input must be a valid file')
end


%    FID = fopen(filename);
%    hapmapTline = fgets(FID);         %read header line
%	word = gettext(Tline);
%	if strcmp(word,'rs#')==1
%	  Tline = fgets(FID);
%	end
%    fclose(FID)

if noise,
	disp(['Reading phased haplotype data file ',filename,' ...'])
end

txt = textread(filename,'%s','delimiter','\n','whitespace','','bufsize',8000);

% find first empty string in cell array, which occurs after the first
% consensus line
mt = find(cellfun('isempty',txt));
% eliminate empty lines
txt(mt) = [];

[yes1,idx1]=ismember('snps:', txt);
[yes2,idx2]=ismember('phased_haplotypes:', txt);

if ~(yes1 && yes2)
      warning('File format error');
      hapldata=[];
      markinfo=[];
      return;
else
      makinftxt=txt(idx1+1:idx2-1);
      haplottxt=txt(idx2+1:end);
end


pos=zeros(1,length(makinftxt));
for (k=1:length(makinftxt)),
      makinf=makinftxt{k};
      makinf=makinf(strfind(makinf,'rs'):end);
      x=strfind(makinf,':');
      rsid{k}=makinf(1:x-1);
      pos(k)=str2double(makinf(x+1:end));
end
markinfo.rsid=rsid;
markinfo.pos=pos;

%hapldata
S='';
pos=zeros(1,length(haplottxt));
for (k=1:length(haplottxt)),
      haplot=haplottxt{k};
      x=strfind(haplot,':');
      Ss=haplot(x+2:end);
      S=strvcat(S,Ss);
end

%hapldata = encodeseq(S);
[hapldata] = i_numeric(S);

markinfo.maf=snp_maf(hapldata,1);


if (size(hapldata,2)~=length(makinftxt))
	error('File format error.');
end





%function [p]=i_maf(hapldata)
%[n,m]=size(hapldata);
%p=zeros(m,1);
%for (k=1:m),
%      x=hapldata(:,k);
%      x(find(sum((x==5),2)>0),:)=[]; % remove
%      [a,b,c]=unique(x(:));
%      y=sum(c==1)/length(x(:));
%      y=min(y,1-y);
%      p(k)=y;
%end


function [genodata] = i_numeric(genodata)
	genodata=double(genodata);
	genodata(find(genodata==double('A')))=1;
	genodata(find(genodata==double('C')))=2;
	genodata(find(genodata==double('G')))=3;
	genodata(find(genodata==double('T')))=4;
	genodata(find(genodata==double('N')))=5;
