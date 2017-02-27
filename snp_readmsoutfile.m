function [hapldata,markinfo]=snp_readmsoutfile(filename,idx)
% [hapldata,markinfo] = snp_readmsoutfile(filename)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<2, idx=1; end
if nargin < 1
    [filename, pathname] = uigetfile( ...
       {'*.txt;*.out', 'MS Output Files (*.txt, *.out)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a MS output format file');
	if ~(filename), hapldata=[]; markinfo=[]; return; end
	filename=[pathname,filename];
end
if isstruct(filename)
    OUT=filename;
else
    [OUT]=readmsoutput(filename);
end

hapldata=OUT.gametes{idx};
markinfo.daf=sum(hapldata)./size(hapldata,1);
%hapldata=hapldata+1;
a=1:length(OUT.positions{idx});
C=textscan(num2str(a),'%s');
markinfo.rsid=C{1};
markinfo.pos=OUT.positions{idx};

