function [ihs, hmarkinfo] = snp_ihsscatter(haplodata, hmarkinfo, showwaitbar)
%SNP_IHS - Integrated haplotype score |iHS| Scatter Plot.
%
%  Syntax: [ihs]=snp_ihsscatter(haplodata,hmarkinfo)

% if nargin<3, regsize=1.0; end
% if nargin<2, popid='CEU'; end
% if nargin<1, rsid='rs7284767'; end
% if strcmp(popid,'HCB'), popid='CHB'; end
% popid=upper(popid);
% if ~(ismember(popid,{'CEU','CHB','JPT','YRI','JPT+CHB'})),
%   error('Not available population code.');
% end

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin < 3
    showwaitbar = false;
end

pmaf = hmarkinfo.maf;
cutoff = 0.015;
idx = find(pmaf >= cutoff);

if (isempty(idx))
    error('No valid SNPs (MAF>=0.015)');
end

if (length(idx) < 3)
    error('Too few valid SNPs (MAF>=0.015)');
end


hmarkinfo.pos = hmarkinfo.pos(idx);
hmarkinfo.maf = hmarkinfo.maf(idx);
hmarkinfo.rsid = hmarkinfo.rsid(idx);

pos = hmarkinfo.pos;
nx = length(pos);

haplodata = haplodata(:, idx);

n = size(haplodata, 2);
ihs = zeros(1, nx);

if showwaitbar, h = waitbar(0, 'Please wait...'); end
for (k = 1:n),
    if showwaitbar, waitbar(k/n, h); end
    ehh = snp_ehh(haplodata, k);
    if (size(ehh, 1) == 2),
        ihh1 = trapz(ehh(1, :));
        ihh2 = trapz(ehh(2, :));
        if ihh2 > 0
            ihs(k) = abs(log(ihh1/ihh2));
        else
            ihs(k) = nan;
        end
    end
end
if showwaitbar, close(h); end

%{
disp('iHS scores were standardized empirically with the CEU HapMap data (mean = 0.7970; s.d.= 0.6040)');
figure;
plot(pos,ihs,'.');
ylabel('|iHS|')
xlabel('Genomic Postion')
title('|iHS| Scatter Plot')
hline(2.0)
%}

%iHS scores can be standardized using estimates of the mean and s.d.
%obtained via coalescent simulation under a variety of demographic models.
%These simulations were tailored to match the frequency spectrum, SNP
%density and recombination profile of the observed data.
%Alternative demographic models included either exponential growth or a
%bottleneck (which varied in onset, severity, duration and population size
%recovery after the bottleneck).

if (nargout < 1)
    i_dispheader('SNP iHS vs. position')
    for (k = 1:nx),
        fprintf('%s\t%d\t%f\n', hmarkinfo.rsid{k}, hmarkinfo.pos(k), ihs(k));
    end
    i_dispfooter
end
