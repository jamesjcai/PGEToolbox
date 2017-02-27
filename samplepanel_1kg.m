function [idvid,popid,conid]=samplepanel_1kg(uselocal)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<1
	uselocal=true;
end
if uselocal
	load('samplepanel_1kg');
else

	%sample_panel_file='ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel';
	sample_panel_file='ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/integrated_call_sets/integrated_call_samples.20101123.ALL.panel';
	a=urlread(sample_panel_file);
	content=textscan(a,'%s','delimiter','\n');

	s=content{1};
	n=length(s);

	if n~=1092, error('xxx'); end

	idvid=cell(n,1); idvid(:)={'NN00000'};
	popid=cell(n,1); popid(:)={'XXX'};
	conid=cell(n,1); conid(:)={'XXX'};

	for k=1:n
	    x=textscan(s{k},'%s');
	    idvid{k}=x{1}{1};
	    popid{k}=x{1}{2};
	    conid{k}=x{1}{3};
	end
end

