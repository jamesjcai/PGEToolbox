function [markinfo1] = snp_mergemarkinfo(markinfo1, markinfo2)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

markinfo1.allele = [markinfo1.allele; markinfo2.allele];
markinfo1.strand = [markinfo1.strand; markinfo2.strand];
markinfo1.rsid = [markinfo1.rsid; markinfo2.rsid];
markinfo1.chr = [markinfo1.chr; markinfo2.chr];
markinfo1.pos = [markinfo1.pos; markinfo2.pos];
