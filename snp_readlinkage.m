function [genodata,markinfo] = snp_readlinkage(filename,varargin)
%SNP_READLINKAGE - Reads linkage pedigree data
% REF: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped
% http://www.broadinstitute.org/science/programs/medical-and-population-gen
% etics/haploview/input-file-formats-0

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-04-07 11:23:42 -0500 (Sun, 07 Apr 2013) $
% $LastChangedRevision: 523 $
% $LastChangedBy: jcai $

if nargin < 1 || isempty(filename)
%     if ispref('PGEToolbox','lastworkingdir')
% 	try
%         cd(getpref('PGEToolbox','lastworkingdir',pwd));
% 	catch ME
% 	rmpref('PGEToolbox','lastworkingdir');	
% 	end
%     end
    [fname, pathname] = uigetfile( ...
       {'*.ped;*.pedigree', 'Linkage Format Files (*.ped, *.pedigree)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a linkage format file');
	if ~(fname), genodata=[]; markinfo=[]; return; end
	%filename=[pathname,filename];
    filename=fullfile(pathname,fname);
end

pnames = {'Delimiter','MissingGenotype','UseACGT','Noise'};
dflts  = {'\t','0','true','true'};

if nargin==1 || nargin<1    
    answer=i_getParamters(dflts);
    answer{3}=strcmpi('true',answer{3});
    answer{4}=strcmpi('true',answer{4});
    dflts=answer;   
    %dflts  = {' '         '0'               true       true};
end
[delimiter,missinggenotype,useacgt,noise] = parseArgs(pnames,dflts,varargin{:});

%if nargin<2, noise=true; end
if ~exist(filename,'file')
    error('Input must be a valid file')
else
    %setpref('PGEToolbox','lastworkingdir',pathname)
end

if noise, fprintf('Reading linkage file %s...\n',filename); end

%{

genodata=[];
block_size = 10000;
while ~feof(fid)
   blockcontent = textscan(fid, '%s', block_size, 'CommentStyle','#',...
                        'Delimiter','\n','BufSize',16777215);
   txtlines=blockcontent{1};
   for k=1:length(txtlines)
        tline=txtlines{k};
        c=regexp(tline,delimiter,'split');
        marker=[c{7:end}];
        if useacgt, marker=i_nt2int(marker); end
        genodata=[genodata; marker];
   end
end

%}

%%

fid = fopen(filename,'r');
tline=fgetl(fid);
if ~ischar(tline), error('xxx'); end
c=regexp(tline,delimiter,'split');
n=length(c(7:end));
m=0;
while ischar(tline)
   if ~isempty(strtrim(tline))
        m = m + 1;
   end
   tline = fgetl(fid);
end

if (mod(n,2)>0)
    fclose(fid);
    error(['Genotype file appears to have an odd number of lines. ',...
           'Each individual is required to have two chromosomes.']);
end



%%
genodata=uint8(zeros(m,n));
familyid=cell(m,1);  familyid(:)={''};
indvid=cell(m,1);    indvid(:)={''};
fatherid=cell(m,1);  fatherid(:)={''};
motherid=cell(m,1);  motherid(:)={''};
sexid=cell(m,1);     sexid(:)={''};
phenotype=cell(m,1); phenotype(:)={''};

frewind(fid);

y=0;
tline=fgetl(fid);

%%
while ischar(tline)
   c=regexp(tline,delimiter,'split');
   y=y+1;
   
   if useacgt
       g=[c{7:end}];
       geno=i_nt2int(g); 
   else
       geno=uint8(cellfun(@str2double,c(7:end)));
   end
   
   genodata(y,:)=geno;
   tline = fgetl(fid);
end
fclose(fid);


markinfo.familyid=familyid;
markinfo.indvid=indvid;
markinfo.fatherid=fatherid;
markinfo.motherid=motherid;
markinfo.sexid=sexid;
markinfo.phenotype=phenotype;


[~,a,b]=snp_maf(genodata);
for k=1:length(a)
    markinfo.allele{k}=sprintf('%d/%d',a(k),b(k));
    markinfo.rsid{k}=sprintf('Marker_%d',k);
    markinfo.strand{k}='+';
    markinfo.chrid{k}='1';
end
markinfo.pos=1:length(a);
markinfo.popid='Unknown';





function  [seq]=i_nt2int(nt)
    nt=lower(nt);
    nt=uint8(nt)+1-uint8('a');
    nt(nt>20|nt<1)=21;
    map=uint8([1 0 2 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 4 5]);
    seq=map(nt);
    

   
function answer=i_getParamters(defaultanswer)

%'Delimiter' 'MissingGenotype' 'UseACGT'
if nargin<1
    defaultanswer={'\t','0','true','true'};
end
prompt={'Delimiter:',...
        'Missing Genotype:',...
        'Use ACGT:',...
        'Noise:'};
       
name='Input for function';
numlines=1;
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';

answer=inputdlg(prompt,name,numlines,defaultanswer,options);      
