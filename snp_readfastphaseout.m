function [hapldata,markinfo] = snp_readfastphaseout(filename,noise)
%snp_readphaseout - Reads fastPHASE output file
%[hapldata,markinfo] = snp_readfastphaseout(filename,noise)

% Population Genetics & Evolution Toolbox, (C) 2009
%% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-02-28 13:55:55 -0600 (Thu, 28 Feb 2013) $
% $LastChangedRevision: 462 $
% $LastChangedBy: jcai $

if nargin<1
    [filename, pathname] = uigetfile( ...
       {'*.phased;*.out', 'fastPHASE Output Files (*.phased, *.out)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a fastPHASE output file');
	if ~(filename), hapldata=[]; markinfo=[]; return; end
	filename=[pathname,filename];
end
if nargin < 2, noise=0; end

hapldata=[];
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
	disp(['Reading fastPHASE output file ',filename,' ...'])
end

%filename='Y:/biotools/phase/out/edaBC.out';
txt = textread(filename,'%s','delimiter','\n','whitespace','');

% find first empty string in cell array, which occurs after the first
% consensus line
mt = find(cellfun('isempty',txt));
% eliminate empty lines
txt(mt) = [];

[y1,idx1]=ismember('BEGIN GENOTYPES',txt);
[y2,idx2]=ismember('END GENOTYPES',txt);

if ~(y1 && y2 )
    error('File format error');
else
    txtpair=txt(idx1+1:idx2-1);
end

n=length(txtpair);
markinfo.sampleid=txtpair(1:3:n);

txtpair(1:3:n)=[];

hapldata=[];
for k=1:length(txtpair)
    %hapldata=[hapldata; str2double(txtpair{k})];
    hapldata=[hapldata; i_nt2int(txtpair{k})];
end
markinfo.pos=1:size(hapldata,2);


function  [seq]=i_nt2int(nt)
    nt=sscanf(nt,'%s');
    %nt
    %pause
    nt=lower(nt);
    nt=uint8(nt)+1-uint8('a');
    nt(nt>20|nt<1)=21;
    map=uint8([1 0 2 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 4 5]);
    seq=map(nt);
