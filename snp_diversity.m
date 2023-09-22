function [h, hv, h_unbiased, hv_unbiased] = snp_diversity(geno)
%SNP_DIVERSITY - SNP diversity
%  Syntax: [h]=snp_diversity(geno)
%
%The estimator for diversity is (1-pi^2+qi^2)[n/(n – 1)],
%
%Same as: snp_heterozygosity, but gives unbiased esitmation too
%
%ref: doi:10.1046/j.1365-294X.2002.01491.x

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2014-03-26 09:30:26 -0500 (Wed, 26 Mar 2014) $
% $LastChangedRevision: 758 $
% $LastChangedBy: jcai $

p = snp_maf(geno);
n = snp_samplen(geno);

q = 1 - p;
Px = 1 - (p.^2 + q.^2);
Nx = n / (n - 1);

%For an unbiased estimate we corrected h by n/(n-1),
%where n is the number of sequenced individuals.
hv = (n / (n - 1)) .* Px;
h = mean(hv);


%Estimates gene diversity per locus and sample using an unbiased estimator
%(see Nei, 1987, eq 7.39 p 164).
if nargout > 2 || nargout < 1
    hok = snp_obshet(geno);
    hv_unbiased = Nx .* (Px - hok ./ (2 * n));
    h_unbiased = mean(hv_unbiased);
end

if (nargout < 1)
    i_dispheader('SNP Diversity')
    fprintf('Mean SNP Diversity (unbiased) = %f\n', h);
    fprintf('Mean SNP Diversity (Nei, 1987, eq 7.39 p 164)= %f\n', h_unbiased);
    i_dispfooter
end


%1-hok./h_unbia
%[geno] = snp_12geno(geno);
%geno2=[];
%n=snp_marklen(geno);
%for k=1:n
%    gen=geno(:,2*k-1:2*k);
%    geno2=[geno2;gen];
%end
%1-abs(mean(hok)./mean(h_unbia))
