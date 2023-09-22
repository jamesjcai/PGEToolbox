function [hapldata, markinfo] = snp_readphaseout(filename, noise)
%snp_readphaseout - Reads PHASE output file
%[hapldata,markinfo] = snp_readphaseout(filename,noise)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2014-03-26 09:30:26 -0500 (Wed, 26 Mar 2014) $
% $LastChangedRevision: 758 $
% $LastChangedBy: jcai $

if nargin < 1
    [filename, pathname] = uigetfile( ...
        {'*.phased;*.out', 'PHASE Output Files (*.phased, *.out)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Pick a PHASE output file');
    if ~(filename), hapldata = [];
        markinfo = [];
        return;
    end
    filename = [pathname, filename];
end
if nargin < 2, noise = 0; end

hapldata = [];
%markinfo=[];

if ~(exist(filename, 'file') || exist(fullfile(cd, filename), 'file')),
    error('InvalidInput', 'Input must be a valid file')
end


%    FID = fopen(filename);
%    hapmapTline = fgets(FID);         %read header line
%	word = gettext(Tline);
%	if strcmp(word,'rs#')==1
%	  Tline = fgets(FID);
%	end
%    fclose(FID)

if noise,
    %disp('Reading PHASE output file ',filename,' ...')
end

%filename='Y:/biotools/phase/out/edaBC.out';
txt = textread(filename, '%s', 'delimiter', '\n', 'whitespace', '', ...
    'bufsize', 409600); % thanks for Jinchuan Xing

% find first empty string in cell array, which occurs after the first
% consensus line
mt = find(cellfun('isempty', txt));
% eliminate empty lines
txt(mt) = [];

% Positions of loci
idx = find(cellfun(@isempty, strfind(txt, 'Positions of loci:')) == 0);
a = txt{idx};
a = a(strfind(a, ':')+1:end);
b = textscan(a, '%d');
markinfo.pos = double(b{1}');
for k = 1:length(markinfo.pos)
    markinfo.rsid{k} = sprintf('Mark%d', k);
end


[y1, idx1] = ismember('BEGIN LIST_SUMMARY', txt);
[y2, idx2] = ismember('END LIST_SUMMARY', txt);
[y3, idx3] = ismember('BEGIN BESTPAIRS_SUMMARY', txt);
[y4, idx4] = ismember('END BESTPAIRS_SUMMARY', txt);

if ~(y1 && y2 && y3 && y4)
    error('File format error');
else
    %      makinftxt=txt(idx1+1:idx2-1);
    %      haplottxt=txt(idx2+1:end);
    txtlist = txt(idx1+1:idx2-1);
    txtpair = txt(idx3+1:idx4-1);

end


pairlist = zeros(length(txtpair), 2);
for (k = 1:length(txtpair)),
    a = txtpair{k};
    c = textscan(a(strfind(a, '(')+1:strfind(a, ')')-1), ...
        '%d%d', 'delimiter', ',');
    pairlist(k, :) = [c{1}, c{2}];
end


haplist = [];
for (k = 1:length(txtlist)),
    a = txtlist{k};
    c = textscan(a, '%s');
    t = i_numeric(c{1}{2});
    haplist = [haplist; t];
end


for k = 1:size(pairlist, 1)
    x = haplist(pairlist(k, 1), :);
    y = haplist(pairlist(k, 2), :);
    hapldata = [hapldata; x];
    hapldata = [hapldata; y];
end


markinfo.maf = snp_maf(hapldata, 1);


%markinfo.maf=snp_maf(hapldata,1);
%if (size(hapldata,2)~=length(makinftxt))
%	error('File format error.');
%end


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
        genodata = double(genodata);
        genodata(find(genodata == double('A'))) = 1;
        genodata(find(genodata == double('C'))) = 2;
        genodata(find(genodata == double('G'))) = 3;
        genodata(find(genodata == double('T'))) = 4;
        genodata(find(genodata == double('N'))) = 5;
        genodata(genodata > 4) = 5;
