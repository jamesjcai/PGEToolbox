function [hapldata,markinfo] = snp_readplemout(filename,noise)
%SNP_READPLEMOUT - Reads PLEM output file
%[hapldata,markinfo] = snp_readplemout(filename,noise)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<1
    [filename, pathname] = uigetfile( ...
       {'*.txt', 'PLEM Output Files (*.txt)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a PLEM output file');
	if ~(filename), hapldata=[]; markinfo=[]; return; end
	filename=[pathname,filename];
end
if nargin < 2, noise=0; end

%hapldata=[];
%markinfo=[];

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
	disp(['Reading PLEM output file ',filename,' ...'])
end

%filename='Y:/biotools/phase/out/edaBC.out';
%txt = textread(filename,'%s','delimiter','\n','whitespace','',...
%      'bufsize',409600);      % thanks for Jinchuan Xing
fid=fopen(filename,'r');
txt=textscan(fid,'%s','delimiter','\n','whitespace','',...
      'bufsize',409600);
txt=txt{1};
fclose(fid);

% find first empty string in cell array, which occurs after the first
% consensus line
%mt=cellfun('isempty',txt);
% eliminate empty lines
%txt(mt)=[];

idx=find(cellfun(@isempty,strfind(txt,'Individual'))==0);
if isempty(idx), error('PLEM output file error.'); end

n=2*numel(idx);
t=textscan(txt{idx(1)+2},'%s');
c=t{1}{1};
m=length(c);

hapldata=false(n,m);
for k=1:numel(idx)
    x=1;
    for kk=idx(k)+2:idx(k)+3
        x=x*-1;
        t=textscan(txt{kk},'%s');
        c=t{1}{1};        
        for i=1:length(c)
            if c(i)=='1'
                if x>0
                    hapldata(k*2,i)=true;
                else
                    hapldata(k*2-1,i)=true;
                end
            end
        end
    end
end


if nargout>1
    try
        % Rows of ID
        idx=find(cellfun(@isempty,strfind(txt,'ID'))==0);
        if isempty(idx)||~isscalar(idx), error('PLEM output file error.'); end
        f=[];
        Gx=[];
        for k=idx+1:length(txt)
            t=textscan(txt{k},'%d%f%f%s');
            f=[f;t{2}];
            c=t{4}{1};    
            g=false(1,length(c));
            for i=1:length(c)
                if c(i)=='1'
                    g(i)=true;
                end
            end
            Gx=[Gx;g]; 
        end
        markinfo.haptype=logical(Gx);
        markinfo.hapfreq=f;
    catch ME
        markinfo=[];
    end
end




