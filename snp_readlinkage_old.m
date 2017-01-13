function [genodata,markinfo] = snp_readlinkage_old(filename,varargin)
%SNP_READLINKAGE - Reads linkage pedigree data
% REF: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped
% http://www.broadinstitute.org/science/programs/medical-and-population-gen
% etics/haploview/input-file-formats-0

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin < 1
    [filename, pathname] = uigetfile( ...
       {'*.pedigree;*.ped', 'Linkage Format Files (*.pedigree, *.ped)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a linkage format file');
	if ~(filename), genodata=[]; markinfo=[]; return; end
	filename=[pathname,filename];
end


noise=0;
if nargin < 2, noise=0; end

if ~(exist(filename,'file') || exist(fullfile(cd,filename),'file')),
    error('InvalidInput','Input must be a valid file')
end

if noise,
%file = fopen(filename, 'r');
disp(['Reading "Linkage" file ',filename,' ...'])
end

%while 1
%   [x,nr] = fscanf(file,'%s\n');
%   if nr == 0 break; end;
%   disp(x)
%end
 %blasttext = fread(file,'*char')';
%fclose(file);
%res=dlmread(filename,'\t')
%res=xlsread(filename)

[pedname,individualid,fatherid,motherid,sex,status,marker]=...
    textread(filename,...
	'%s%s%s%s%s%s%[ 01234\t]%*[^\n]',...
	'delimiter','\t','commentstyle','shell','bufsize',409600);

%markinfo= struct([]);
markinfo.pedigreename=pedname;
markinfo.individualid=individualid;
markinfo.fatherid=fatherid;
markinfo.motherid=motherid;
markinfo.sex=sex;
markinfo.affectionstatus=status;


caseno=length(marker);
genodata=[];
for k=1:caseno
      genodata = [genodata; sscanf(marker{k},'%d')'];
end

m=size(genodata,2);
if m<3
    error('Genotype file appears to have fewer than 3 columns.')
end
if (mod(m,2)>0),
    error(['Genotype file appears to have an odd number of lines. ',...
           'Each individual is required to have two chromosomes.']);
end
genodata(genodata<1 | genodata>4)=5;



[~,a,b]=snp_maf(genodata);
for k=1:length(a)
    markinfo.allele{k}=sprintf('%d/%d',a(k),b(k));
    markinfo.rsid{k}=sprintf('Marker_%d',k);
    markinfo.strand{k}='+';
    markinfo.chrid{k}='1';
end
markinfo.pos=1:length(a);
markinfo.popid='Unknown';


% function [p] = predhet(genodata,fatherid,motherid)
% [n,m]=size(genodata); %assert(m==2)
%
% count=zeros(1,4);
% for (k=1:n),
%       allele1=genodata(k,1);
%       allele2=genodata(k,2);
%       if (allele1~=0 &&allele2~=0)
% 	if (fatherid(k)==0 && motherid(k)==0)    % no Ancestor
% 	count(allele1)=count(allele1)+1;
% 	count(allele2)=count(allele2)+1;
% 	end
%       end
% end
%
% count(find(count==0))=[];
% p=1-sum(count.^2)/(sum(count)^2);

