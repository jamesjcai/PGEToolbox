function [genodata,markinfo] = snp_readvcf(filename,noise)
%SNP_READVCF - Reads VCF (Variant Call Format) version 4.1
%[genodata,markinfo] = snp_readhapmap(filename,noise)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-12-27 00:25:00 -0600 (Fri, 27 Dec 2013) $
% $LastChangedRevision: 755 $
% $LastChangedBy: jcai $

if nargin < 1
    [filename, pathname] = uigetfile( ...
       {'*.vcf', 'Variant Call Format Files (*.vcf)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a Variant Call Format (VCF) file');
	if ~(filename), genodata=[]; markinfo=[]; return; end
	filename=[pathname,filename];
end

if nargin < 2, noise=false; end

genodata=[];
markinfo=struct();

if ~(exist(filename,'file') || exist(fullfile(cd,filename),'file')),
    error('SNP_READVCF:argChk','Input must be a valid file')
end


    fid = fopen(filename,'r');
%    hapmapTline = fgets(FID);         %read header line
%	word = gettext(Tline);
%	if strcmp(word,'rs#')==1
%	  Tline = fgets(FID);
%	end
%    fclose(FID)

if noise, fprintf('Reading VCF file %s...\n',filename); end
ignc=0;
iclc=0;
G=[];

block_size = 10000;

try
   % while ~feof(fid)
       s = textscan(fid, '%s', block_size, 'CommentStyle','#',...
                            'Delimiter','\n','BufSize',16777215);
       s=s{1};
       
       if isempty(s)
           fclose(fid);
           return;
       else
           for k=1:length(s)
                r=regexp(s{k},'\t','split');
                maja=upper(r{4});
                mina=upper(r{5});
                if ismember(maja,{'A','C','G','T'})&&ismember(mina,{'A','C','G','T'})
                    iclc=iclc+1;
                    chrid{iclc}=r{1};
                    pos(iclc)=uint32(str2double(r{2}));
                    rsid{iclc}=r{3};
                    refalle(iclc)=i_nt2int(maja);
                    altalle(iclc)=i_nt2int(mina);            
                    g=i_process_geno(r(10:end),refalle(iclc),altalle(iclc));
                    %size(G)
                    %size(g)

                    G=[G; g];
                else
                    ignc=ignc+1;
                end
           end
       end
    % end
catch ME
    disp(ME.message);
    fclose(fid);
    return;
end
fclose(fid);
genodata=snp_genotranspose(G);
markinfo.rsid=rsid';
markinfo.refalle=refalle';
markinfo.altalle=altalle';
markinfo.pos=pos';
markinfo.chrid=chrid';



function g=i_process_geno(x,refalle,altalle)
g=uint8(5*ones(1,length(x)*2));
for k=1:length(x)    
    if any(strcmpi(x{k}(1:3),{'0|0','0/0'}))
        g([k*2-1,k*2])=[refalle refalle];
    elseif any(strcmpi(x{k}(1:3),{'1|1','1/1'}))
        g([k*2-1,k*2])=[altalle altalle];
    elseif any(strcmpi(x{k}(1:3),{'1|0','1/0'}))
        g([k*2-1,k*2])=[altalle refalle];
    elseif any(strcmpi(x{k}(1:3),{'0|1','0/1'}))
        g([k*2-1,k*2])=[refalle altalle];
    end
end



function  [seq, map] = i_nt2int(nt)
if isempty(nt)
    seq = uint8([]);
    return
end

unknown = 0;
acgtonly = true;
origsize = size(nt);
nt = nt(:);
nt = nt';

% * A C G T U R Y K M S  W  B  D  H  V N
% 0 1 2 3 4 4 5 6 7 8 9 10 11 12 13 14 15

%           'a  b c  d e f g  h i j k l m  n o p q r s t u  v  w  x y z  -'
map = uint8([1 11 2 12 0 0 3 13 0 0 7 0 8 15 0 0 0 5 9 4 4 14 10 15 6 0 16]);
if acgtonly
    map(map > 4) = 0;
end
if unknown ~= 0
    map(map == 0) = unknown;
end
nt = lower(nt);
seq = map((uint8(nt) + 1) - uint8('a'));
seq = reshape(seq,origsize);

