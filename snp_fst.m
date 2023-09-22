function [fv] = snp_fst(geno1, geno2, methodtag)
%SNP_FST - Fst SNPs from two populations
%
%    [fst]=snp_fst(geno1,geno2,'Weir')
%    [fst]=snp_fst(geno1,geno2,'Wright')
%    [fst]=snp_fst(geno1,geno2,'Hughes')
%
% See also: SNP_FST3

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-30 18:59:03 -0600 (Wed, 30 Jan 2013) $
% $LastChangedRevision: 380 $
% $LastChangedBy: jcai $

if nargin < 3, methodtag = 'weir'; end

%for x=1:4
%    a=[2*x-1,2*x];
%    snp_fst(geno1(:,a),geno2(:,a))
%end
%snp_fst(geno1(:,x),geno2(:,x))
%0.4869
% my     hiploid  diploid
%0.1980  0.201336   0.202195
%0.5318  0.534956   0.535098
%0.0569  0.062531   0.062883
%0.4832  0.486359   0.486375

if nargin < 1
    disp('Example data')
    load('example_data/twogeno_Fst');
end

[m1] = snp_marklen(geno1);
%[n1,m1]=snp_samplen(geno1);
%[n2,m2]=snp_samplen(geno2);
%f=i_fsthiploid(geno1,geno2);
%if nargout>1
fv = zeros(1, m1);
for k = 1:m1
    genoa = geno1(:, 2*k-1:2*k);
    genob = geno2(:, 2*k-1:2*k);
    switch (lower(methodtag))
        case 'weir'
            fv(k) = i_fst_weir_hiploid(genoa, genob);
        case 'wright'
            fv(k) = i_fst_wright(genoa, genob);
        case 'hughes'
            fv(k) = i_fst_hughes05(genoa, genob);
        otherwise
            error('Wrong METHODTAG');
    end
end


    function [f] = i_fst_wright(geno1, geno2)
        %p1_i = allele freq in pop 1
        %p2_i = allele freq in pop 2
        %i = number of sites

        %REF: http://www.ajhg.org/AJHG/fulltext/S0002-9297(07)63770-7
        %but see: http://intl.genome.org/cgi/content/full/15/11/1496

        [p1, p2, ~, ~, pbar] = i_subpopallelefreq(geno1, geno2);
        q1 = 1 - p1;
        q2 = 1 - p2;
        qbar = 1 - pbar;
        if pbar == 0 || qbar == 0
            f = nan;
        else
            f = 1 - mean([2 * p1 .* q1, 2 * p2 .* q2]) ./ (2 * pbar .* qbar);
        end


        %function [f]=i_st_wahlund(geno1,geno2)
        %http://www.nature.com/ng/journal/v39/n12/pdf/ng.2007.13.pdf
        %Fst=var(p)/p*(1-p*) was used to measure allele frequency differences
        %between populations, where var(p) represents the variance of the
        %frequencies of an allele from a biallelic SNP, and p* represents the
        %average frequency of the allele among the populations under consideration.


        %function [f]=i_fst_hudson92(geno1,geno2)
        %Kst* Hudson et al 1992 an analog of Fst.
        % A test for detecting geographic subdivision MBE


            function [f] = i_fst_hughes05(geno1, geno2)
                %Ref: Hughes et al. 2005
                %p1 and p2 are the frequencies of the first allele in each of two subpopulations
                %http://www.genetics.org/cgi/reprint/170/3/1181
                [p1, p2] = i_subpopallelefreq(geno1, geno2);
                q1 = 1 - p1;
                q2 = 1 - p2;
                f = 1 - sum([sqrt(p1.*p2), sqrt(q1.*q2)]);

                    function [f] = i_fst_weir_hiploid(geno1, geno2)
                        [p1, p2, n1, n2] = i_subpopallelefreq(geno1, geno2);
                        if p1 == p2
                            f = 0;
                        else
                            f = fst_weir(n1, n2, p1, p2);
                        end
                        %REF: http://mbe.oxfordjournals.org/cgi/content/full/23/9/1697

                            function [p1, p2, n1, n2, pbar] = i_subpopallelefreq(geno1, geno2)

                                %s=2;                   % for SNP pairs, num of subpoulations, s=2
                                n1 = sum(geno1(:) ~= 5); % hiploid so n1 = n1*2;
                                n2 = sum(geno2(:) ~= 5);
                                %n=n1+n2;

                                %nc = (1/(s-1))*((n1+n2)-(n1^2+n2^2)/(n1+n2));
                                %nc = (n-(n1^2+n2^2)/n);

                                [p1, alle1] = snp_maf(geno1);
                                [p2, alle2] = snp_maf(geno2);
                                if alle1 ~= alle2
                                    p2 = 1 - p2;
                                end
                                if nargout > 4
                                    [pbar, ~] = snp_maf([geno1; geno2]);
                                end
