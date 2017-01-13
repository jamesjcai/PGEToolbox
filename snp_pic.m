function [pic]=snp_pic(geno)
%SNP_PIC - SNP Polymorphism Information Content
%
%  Syntax: [pic]=snp_pic(geno)
%
% The PIC value is defined as the probability that a given marker genotype
% of an offspring of an affected parent will allow deduction of the
% parental genotype at the marker locus (Botstein et al. 1980)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if any(size(geno)==1), 
    p=geno; % input is the allele frequency
else
    p=snp_maf(geno);
end
q=1-p;

pic=1-(p.^2+q.^2)-2*(p.^2).*(q.^2);

if (nargout<1)
i_dispheader('SNP Polymorphism Information Content')
	fprintf ('PIC = %g\n',pic);
i_dispfooter
end

