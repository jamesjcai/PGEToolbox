function [genodata,markinfo] = snp_readhapmap(filename,noise)
%SNP_READHAPMAP - Reads HapMap genotype data
%[genodata,markinfo] = snp_readhapmap(filename,noise)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin < 1
    [filename, pathname] = uigetfile( ...
       {'*.hapmap;*.hmp', 'HapMap Files (*.hapmap, *.hmp)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a HapMap (SNP genotype data) format file');
	if ~(filename), genodata=[]; markinfo=[]; return; end
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
	disp(['Reading "HapMap" file ',filename,' ...'])
end

[rsid,alleles,chrom,pos,strand,assemblyid,center,...
protLSID,assayLSID,panelLSID,QCcode,marker] = textread(filename,...
	'%s%s%s%d%s%s%s%s%s%s%s%[ATCGN ]%*[^\n]','delimiter',' ',...
	'commentstyle','shell','headerlines',3);

%markerinfo= struct('rsid',rsid,'alleles',alleles,'chrom',chrom,'pos',pos,'strand',strand,'assemblyid',assemblyid,'center',center,...
%'protLSID',protLSID,'assayLSID',assayLSID,'panelLSID',panelLSID,'QCcode',QCcode);

n=length(marker);
if n>0
    popid=i_getpopid(panelLSID);
else
    popid=[];
end
genodata=[];
for (k=1:n),
%      marker{k}
%      pause
%      sscanf(marker{k},'%s')
%      pause
      genodata = [genodata; sscanf(marker{k},'%s')];
end
genodata = i_hapmap_transpose(genodata);
genodata = i_numeric(genodata);


markinfo.allele=alleles;
markinfo.strand=strand;
markinfo.rsid=rsid;
markinfo.chr=chrom;
markinfo.pos=pos;
markinfo.popid=popid;



%indvinfo.fatherid=zeros(1,size(genodata,1));
%indvinfo.motherid=zeros(1,size(genodata,1));




function [geno_out] = i_hapmap_transpose(geno_in)
	%geno_in=['GCGGGC';'AAATAT'];
	[n,m]=size(geno_in);
	%geno_out=char(ones(m/2,n*2)*42);
	geno_out=zeros(m/2,n*2);

	for (j=1:n),
	for (k=1:2:m),
	      x=(k+1)/2;
	      y=j*2-1;
	      geno_out(x,[y y+1]) = geno_in(j,[k k+1]);
	end
	end


function [genodata] = i_numeric(genodata)
	genodata=double(genodata);
	genodata(genodata==double('A'))=1;
	genodata(genodata==double('C'))=2;
	genodata(genodata==double('G'))=3;
	genodata(genodata==double('T'))=4;
	genodata(genodata==double('N'))=5;
    genodata(genodata<1 | genodata>4)=5;


function popid=i_getpopid(panelLSID)
  if isempty(panelLSID)
            popid=[];
  else
            panelLSID=panelLSID{1};
            if strcmp(panelLSID,'Yoruba')
                popid='YRI';
            elseif strcmp(panelLSID,'Yoruba')
                popid='CEU';
            elseif strcmp(panelLSID,'Japanese+Han')
                popid='JPN+CHB';
            elseif strcmp(panelLSID,'Japanese')
                popid='JPT';
            elseif strcmp(panelLSID,'Han')
                popid='CHB';
            else
                popid='Unknown';
            end
%         for k=1:length(panelLSID)
%             if strcmp(panelLSID{k},'Yoruba')
%                 popid{k}='YRI';
%             elseif strcmp(panelLSID{k},'Yoruba')
%                 popid{k}='CEU';
%             elseif strcmp(panelLSID{k},'Japanese+Han')
%                 popid{k}='JPN+CHB';
%             elseif strcmp(panelLSID{k},'Japanese')
%                 popid{k}='JPT';
%             elseif strcmp(panelLSID{k},'Han')
%                 popid{k}='CHB';
%             else
%                 popid{k}='Unkonwn';
%             end
%         end

  end