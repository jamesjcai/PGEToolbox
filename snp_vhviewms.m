function [G] = snp_vhviewms(msout)
%SNP_VHVIEWMS - visual haplotype for MS output
% Syntax: snp_vhviewms(msout)
%
%    blue =  common allele
%    yellow = homozygous genotype for the rare allele
%    gray - missing data (N)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-04-23 23:03:55 -0500 (Tue, 23 Apr 2013) $
% $LastChangedRevision: 530 $
% $LastChangedBy: jcai $


if isempty(msout), return; end
G = [];
[n, m] = size(msout);
msout = msout + 1;
G = cat(2, msout, ones(n, 1)*3);
G = cat(1, G, ones(1, m+1)*3);
G(G > 3) = 3;

pcolor(G);
axis ij;
grid off;

colormap('default');
ax = jet;
colormap(ax([1, 40, 24], :))

axis equal
%axis off
axis tight

return;

x = axis;
t = max(x(2), x(4));
x(2) = t;
x(4) = t;
axis(x);
axis equal
%axis off
axis tight
xlabel('Markers (SNPs)')
ylabel('Samples (Chromosomes)')


i_dispheader('Visual Haplotype (VH) View')
disp('Blue   : Homozygote-common allele ')
disp('Yellow : Homozygote-rare allele')
disp('Cyan   : Undetermined ')
disp(' ')
disp('* Column = Markers (SNPs)')
disp('** Row = Samples (Chromosomes)')
i_dispfooter
%shading flat
