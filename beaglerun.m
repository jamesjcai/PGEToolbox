function [out] = beaglerun(genodata, markeinfo)
%BEAGLERUN - run BEAGLE to make haplotype clusters
%
% [ldinfo]=beaglerun(genodata)
%
%SEE ALSO:

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-04-23 19:44:14 -0500 (Tue, 23 Apr 2013) $
% $LastChangedRevision: 529 $
% $LastChangedBy: jcai $

if nargin < 2
    markeinfo = [];
end
oldpath = pwd;
out = [];
cdpge;
cd 'addins/EMLD';
%[exedir,dlgshown]=pge_getprgmdir(sprintf('%s_prgmdir',mfilename));
%if isempty(exedir)||dlgshown, return; end
%cd(exedir);

filename = 'input.bgl';
try
    disp('Writing BEAGLED input file...')
    snp_writebeagle(genodata, markeinfo, filename);
    disp('Running java EMLD...')
    system(sprintf('java -Xmx500m -jar beagle.jar data=%s out=output', filename));
    % java -Xmx500m -jar beagle.jar unphased=zzz.bgl missing=? out=examplex
catch
    out = [];
    return;
end
disp('Reading BEAGLE output...')
gunzip('output.input.bgl.dag.gz')
[out] = textread('output.input.bgl.dag', '%s', 'headerlines', 1, 'delimiter', '\n');
cd(oldpath);
