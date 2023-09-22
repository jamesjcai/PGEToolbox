function snp_viewmark(markinfo)
%SNP_VIEWMARK - Displays SNP marker information

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

alleles = markinfo.allele;
strands = markinfo.strand;
rsids = markinfo.rsid;
position = markinfo.pos;
chrom = markinfo.chr;

i_dispheader('View Marker Information')
fprintf('rsid\tAllele\tChromosome\tStrand\tPosition\n');
for (k = 1:length(position)),
    fprintf('%10s\t%s\t%s\t%s\t%s\t%d\n', ['Mrk_', num2str(k)], rsids{k}, alleles{k}, chrom{k}, strands{k}, position(k));
end
if isfield(markinfo, 'popid')
    fprintf('==========================================\n');
    fprintf('Population: %s\n', markinfo.popid);
end
i_dispfooter