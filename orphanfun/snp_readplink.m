function [genodata,markinfo] = snp_readplink(filename,noise)
%SNP_READLINKAGE - Reads linkage pedigree data

% Population Genetics & Evolution Toolbox, (C) 2007
% Author: James J. Cai
% Email: jamescai@stanford.edu
% Website: http://bioinformatics.org/pgetoolbox/
% Last revision: 2/23/2007



%The autosomes should be coded 1 through 22. 
%The following other codes can be used to specify other chromosome types:
%
%     X    X chromosome                    -> 23
%     Y    Y chromosome                    -> 24
%     XY   Pseudo-autosomal region of X    -> 25
%     MT   Mitochondrial                   -> 26


%     Family ID
%     Individual ID
%     Paternal ID
%     Maternal ID
%     Sex (1=male; 2=female; other=unknown)
%     Phenotype
     
     
if nargin < 1
    [filename, pathname] = uigetfile( ...
       {'*.pedigree;*.ped', 'PED Files (*.pedigree, *.ped)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a PED file');
	if ~(filename), genodata=[]; markinfo=[]; return; end
	filename=[pathname,filename];
end

if nargin < 2, noise=0; end

if ~(exist(filename,'file') || exist(fullfile(cd,filename),'file')),
    error('InvalidInput','Input must be a valid file')
end

if noise, 
%file = fopen(filename, 'r');
disp(['Reading PED file ',filename,' ...'])
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

try
[pedname,id,fatherid,motherid,sex,status,marker] = textread(filename,...
	'%s%n%n%n%n%n%[ 01234\t]%*[^\n]',...
	'delimiter','\t','commentstyle','shell');
catch
try
    [pedname,id,fatherid,motherid,sex,status,marker] = textread(filename,...
	'%n%n%n%n%n%n%s[^\n]',...
	'delimiter',' ','commentstyle','shell');
catch
    disp('yyy')
end
end

%markinfo= struct([]);
markinfo.familyid = pedname;
markinfo.fndividualid=id;
markinfo.fatherid= fatherid;
markinfo.motherid= motherid;
markinfo.sex=sex;
markinfo.phenotype= status;


caseno=length(marker);
genodata=[];
for (k=1:caseno),
      genodata = [genodata; sscanf(marker{k},'%d')'];
end

[n,m]=size(genodata);
if (m<3),
error('Genotype file appears to have fewer than 3 columns.')
end
if (mod(m,2)>0),
error('Genotype file appears to have an odd number of lines. Each individual is required to have two chromosomes.')
end

markinfo


function [p] = predhet(genodata,fatherid,motherid)
[n,m]=size(genodata); %assert(m==2)

count=zeros(1,4);
for (k=1:n),
      allele1=genodata(k,1);
      allele2=genodata(k,2);
      if (allele1~=0 &&allele2~=0)
	if (fatherid(k)==0 && motherid(k)==0)    % no Ancestor
	count(allele1)=count(allele1)+1;
	count(allele2)=count(allele2)+1;
	end
      end
end


count(find(count==0))=[];
p=1-sum(count.^2)/(sum(count)^2);