function [h,hvar]=snp_haphet(hap)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


[numHap,sizHap,seqHap]=counthaplotype(hap);
p=sizHap./sum(sizHap);
h=1-sum(p.^2)./(1-1/numHap);

if nargout>1
    n=numHap;
    hvar=(2*(n-1)/n^3)*(2*(n-2)*(sum(p.^3)-(sum(p.^2))^2));
%(a slight modification of the standard diploid variance; Nei 1987).
% Nei M Molecular evolutionary genetics. In Columbia University Press 1987
% New York:Columbia University Press as Nash D, Nair S, Mayxay M, Newton
% PN, Guthmann JP, Nosten F, Anderson TJ.
%
%Selection strength and hitchhiking around two anti-malarial resistance
%genes. Proc Biol Sci. 2005 Jun 7;272(1568):1153-61. PubMed PMID: 16024377;
%PubMed Central PMCID: PMC1559806.

end