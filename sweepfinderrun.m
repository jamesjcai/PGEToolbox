function [stats] = sweepfinderrun(geno, mark, gridsize, linestyle)
%SWEEPFINDERRUN - executes SWEEPFINDER
%
%sweepfinderrun(geno,mark)
%[stats]=sweepfinderrun(geno,mark,gridsize)
%
%     stats.position - locations along a chromosome
%     stats.clr      - composite likelihood ratio
%     stats.alpha    - the strength of the sweep
%
%"This program (SWEEPFINDER) reads in SNP data with locations along a chromosome.
%It creates a grid of locations over this region, and at each location, calculates
%the maximum composite likelihood ratio (CLR) statistic comparing the hypothesis
%of a complete selective sweep at the location to the null hypothesis of no sweep.
%The likelihiood function is the "parameteric approach" described in the reference
%(eq 6). It outputs the CLR statistic as well as the parameter alpha (the strength
%of the sweep) to an outfile.  Small values of alpha correspond to strong sweeps."
%-- SweepFinder program (Last updated: 8/14/06) Melissa Hubisz and Rasmus Nielsen
% Refererence: Nielsen et al, Genome Research 2005 Nov; 15(11):1566-75

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-04-23 19:44:14 -0500 (Tue, 23 Apr 2013) $
% $LastChangedRevision: 529 $
% $LastChangedBy: jcai $

if nargin < 4
    linestyle = '-';
end
if nargin < 3
    gridsize = 50;
end

%[exedir,dlgshown]=pge_getprgmdir(sprintf('%s_prgmdir',mfilename));
%if isempty(exedir)||dlgshown, return; end

oldpath = pwd;

cdpge;
cd('addins/sweepfinder');
%cd(exedir);
[status] = snp_writesweepfinder(geno, mark, 'input.tab');
if status ~= 1
    cd(oldpath);
    error('Error writing SweepFinder input file.');
end

cmd = sprintf('SweepFinder -s %d input.tab output.tab', gridsize);
%fprintf('Running: %s\n\n',cmd);
system(cmd);

[vpos, vclr, valpha] = textread('output.tab', '%f%f%f', 'headerlines', 1);
if nargout < 1
    vpos = vpos ./ 1000000;
    %figure;
    %subplot(2,1,1)
    %plot(vpos,vclr,'LineSmoothing','on');
    plot(vpos, vclr, linestyle);
    xlim([min(vpos(:)), max(vpos(:))]);
    %ylabel('CLR (Composite Likelihood Ratio)')
    ylabel('CLR')
    xlabel('Position (Mb)')
else
    stats = struct;
    stats.position = vpos;
    stats.CLR = vclr;
    stats.alpha = valpha;
end
cd(oldpath);
