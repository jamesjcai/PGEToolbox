function snp_sfsmwu()
%SNP_SFSMWU - Mann-Whitney U-test (MWU) applied to the frequency spectrum
%
%REF:
%
%Akashi, H. 1999. Inferring the fitness effects of DNA mutations from
%patterns of polymorphism and divergence: Statistical power to detect
%directional selection under stationarity and free recombination. Genetics
%151: 221–238.
%
%http://www.genome.org/cgi/content/full/15/11/1566
%
%Mann-Whitney U-test (MWU) applied to the frequency spectrum can be used to
%test for an excess or deficiency of low-frequency derived alleles (Akashi
%1999).
%
%MWU should be performed using the folded, and not the unfolded frequency
%spectrum, because use of the unfolded frequency spectrum involves
%additional assumptions regarding the outgroup species and may be sensitive
%to mis-specifications of ancestral states ().

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $
