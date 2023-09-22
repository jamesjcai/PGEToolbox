function [F] = snp_fstmat(query, ishapmap3)
%SNP_FSTMAT - returns mean Fst matrix in query region b/w populations
%
%[F]=snp_fstmat(query,ishapmap3)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin < 1, query = 'EDAR'; end
if nargin < 2, ishapmap3 = true; end
if ishapmap3
    %popset={'ASW','CEU','CHB','CHD','GIH','JPT','LWK','MEX','MKK','TSI','YRI'};
    popset = {'LWK', 'MKK', 'YRI', 'ASW', 'CEU', 'TSI', 'GIH', 'CHB', 'CHD', 'JPT', 'MEX'};
else
    popset = {'CEU', 'CHB', 'JPT', 'YRI'};
end

n = length(popset);
F = ones(n) * -1;
for i = 1:n - 1
    for j = i + 1:n
        pop1 = popset{i};
        pop2 = popset{j};
        [fv] = snp_fstquery(query, pop1, pop2, ishapmap3);
        pause(5)
        F(i, j) = nanmean(fv);
        F(j, i) = F(i, j);
    end
end

if nargout < 1

    ttxt = sprintf('Average F_{ST} of SNPs at %s', query);
    %subplot(1,2,1)

    F(isnan(F)) = 0;
    matrixcircle(triu(F, 1), popset, ttxt);
    hline([4.5, 7.5, 10.5], 'r-')
    vline([4.5, 7.5, 10.5], 'r-')
    %subplot(1,2,2)
    %matrixcolor(F,popset,ttxt);
    % hline([5 8 11],'g-')
    % vline([5 8 11],'g-')


    if ishapmap3
        i_dispheader(sprintf('F_{ST} between HapMap3 populations at %s', query));

        fprintf('LWK (L): Luhya in Webuye, Kenya\n');
        fprintf('MKK (K): Maasai in Kinyawa, Kenya\n');
        fprintf('YRI (Y): Yoruba in Ibadan, Nigeria (West Africa)\n');
        fprintf('ASW (A): African ancestry in Southwest USA\n');
        fprintf('CEU (C): Utah residents with Northern and Western European ancestry from the CEPH collection\n');
        fprintf('TSI (T): Toscans in Italy\n');
        fprintf('GIH (G): Gujarati Indians in Houston, Texas\n');
        fprintf('CHB (H): Han Chinese in Beijing, China\n');
        fprintf('CHD (D): Chinese in Metropolitan Denver, Colorado\n');
        fprintf('JPT (J): Japanese in Tokyo, Japan\n');
        fprintf('MEX (M): Mexican ancestry in Los Angeles, California\n');


    else
        i_dispheader(sprintf('F_{ST} between HapMap populations at %s', query));
        fprintf('CEU (C): Utah residents with Northern and Western European ancestry from the CEPH collection\n');
        fprintf('CHB (H): Han Chinese in Beijing, China\n');
        fprintf('JPT (J): Japanese in Tokyo, Japan\n');
        fprintf('YRI (Y): Yoruba in Ibadan, Nigeria (West Africa)\n');
    end
    i_dispfooter
    %set(gca,'XTickLabel',popset);
end
