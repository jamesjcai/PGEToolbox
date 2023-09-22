function [d, thetw, pval] = snp_tajima89d(geno)
%SNP_TAJIMA89D - Calculates Tajima's D based on SNP frequencies
% Syntax: [d,theta,pval] = snp_tajima89d(geno)

% smpln is the number of chromosomes resequenced.  default is 30
% Calculate Tajima's D based on SNP frequencies
% Calculate pi using frequencies of alternative alleles

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


[smpln] = snp_samplen(geno);

[thepi] = snp_thetapi(geno);
[Sn] = snp_segsites(geno);


[d, thetw, pval] = tajima89d(smpln, Sn, thepi);
% thanks for Jinchuan Xing for bug report

if nargout < 1
    i_dispheader('Tajima''s Neutrality Test')
    disp('Mode: SNP');
    fprintf('\n');
    fprintf('No. of Segregating sites (S): %d\n', Sn);
    fprintf('Nucleotide diversity (per sequence), Pi: %f\n', thepi);
    fprintf('\n');
    %fprintf (['Diff = %f, s.e. = %f\n'], Diff, DiffSE);
    fprintf('Tajima''s D: %f\n', d);
    [D] = tajima89d_simu(smpln, 1000, thetw, 0);
    p = sum(d < D) ./ 1000;
    fprintf('Statistical significance:\n P = %f%s\n', p, sigtag(p));
    i_dispfooter
end
