function [D_raw, D_prime, R2, R, D_dist, siteidx, pvalue, done] = linkdisequ(seq, squred, warninglarge)
%LINKDISEQU - Linkage disequilibrium from sequences
%
% Syntax: [D_raw,D_prime,R2,R,D_dist,siteidx] = linkdisequ(seq,squred)
%
% squred - 1 = outputs in squre matrices; 0 (defult) = in vector

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (nargin < 3), warninglarge = 0; end
if (nargin < 2), squred = 0; end
if (isstruct(seq)), seq = seq.seq; end

[n, m] = size(seq);
siteidx = zeros(1, m);
n_bialle = 0;

for k = 1:m
    site = seq(:, k);
    if (length(unique(site)) == 2),
        n_bialle = n_bialle + 1;
        siteidx(n_bialle) = k;
    end
end

if (n_bialle < 2),
    error('No pairwise comparisons.')
end

siteidx(siteidx == 0) = [];
countpair = nchoosek(n_bialle, 2);

D_raw = zeros(n_bialle);
D_prime = zeros(n_bialle);
R2 = zeros(n_bialle);
R = zeros(n_bialle);
D_dist = zeros(n_bialle);
pvalue = zeros(n_bialle);
done = 0;

if (warninglarge && n_bialle > 200)
    answer = questdlg('Do you want to continue?', ...
        'Large number of pairwise comparisons!', ...
        'Continue', 'Cancel', 'Continue');
    switch (lower(answer))
        case 'cancel'
            %helpdlg('Action cancelled.')
            return;
    end
end

if (nargout < 1), fprintf(['\n   site1    site2     Dist          D       D''      R^2        R\n']); end

for i = 1:n_bialle - 1
    for j = i + 1:n_bialle
        site1 = seq(:, siteidx(i));
        site2 = seq(:, siteidx(j));
        [D_raw(i, j), D_prime(i, j), R2(i, j), R(i, j), pvalue(i, j)] = i_ld(site1, site2);
        if ~(nargout > 0 && nargout < 4)
            D_dist(i, j) = abs(siteidx(i)-siteidx(j));
            if (nargout < 1),
                fprintf(['%8d %8d %8d %10.4f %8.4f %8.4f %8.4f\n'], ...
                    siteidx(i), siteidx(j), D_dist(i, j), D_raw(i, j), ...
                    D_prime(i, j), R2(i, j), R(i, j));
            end
        end
    end
end

if (nargout < 1), fprintf(['   site1    site2     Dist          D       D''      R^2        R\n\n\n']); end

if (squred ~= 1),
    D_raw = i_mat2vec(D_raw, n_bialle);
    D_prime = i_mat2vec(D_prime, n_bialle);
    R2 = i_mat2vec(R2, n_bialle);
    R = i_mat2vec(R, n_bialle);
    D_dist = i_mat2vec(D_dist, n_bialle);
    pvalue = i_mat2vec(pvalue, n_bialle);
end


if (nargout < 1),
    i_dispheader('Linkage Disequilibrium')
    fprintf('Number of polymorphic sites analyzed: %d\n', n_bialle);
    fprintf('Number of pairwise comparisons: %d\n', countpair);
    disp(' ')

    % ZnS statistic (Kelly 1997, equation 3). ZnS is the average of R^2 (Hill and
    % Robertson 1968) over all pairwise comparisons.
    ZnS = sum(i_mat2vec(R2, n_bialle)) / countpair;
    fprintf('Value of ZnS (Kelly 1997): %f\n', ZnS);

    %fprintf('Value of Za (Rozas et al. 2001): 0.3202\n');  %Za and ZZ statistics (Rozas et al. 2001). Za is the average of R^2 (Hill and Robertson 1968) over all pairwise comparisons between adjacent polymorphic sites; ZZ = Za - ZnS. ZZ statistic could be used for detecting intragenic recombination.
    %fprintf('Value of ZZ (Rozas et al. 2001): 0.1165\n');
    %fprintf('Value of Wall''s B: 0.2061');
    %fprintf('Value of Wall''s Q: 0.2727');

    disp(' ')
    disp(' ===== Regression Equation:  Y = a + bX  (X measured in kb) =====')
    disp(' ')
    [ps] = polyfit(D_dist, abs(D_raw), 1);
    fprintf([' |D| values:  Y = %f + %fX   (%d points)\n'], ps(2), 1000*ps(1), countpair);
    [ps] = polyfit(D_dist, abs(D_prime), 1);
    fprintf([' |D''| values:  Y = %f + %fX   (%d points)\n'], ps(2), 1000*ps(1), countpair);
    [ps] = polyfit(D_dist, R2, 1);
    fprintf([' r^2 values:  Y = %f + %fX   (%d points)\n'], ps(2), 1000*ps(1), countpair);
    disp(' ')
    %figure;
    %plot(D_dist, D_prime, '.')
    %xcount=length(siteidx);
    %counter=0;
    %for (i=1:xcount-1),
    %for (j=i+1:xcount),
    %	counter=counter+1;
    %	fprintf(['%d %d %f'], siteidx(i), siteidx(j), D_raw(i,j));
    %end
    %end
    i_dispfooter
end
done = 1;


    function [V] = i_mat2vec(M, n_bialle)
        if (sum(size(M) == 1) > 0), V = M;
            return;
        end
        xM = [];
        for k = 1:n_bialle
            xM = [xM, M(k, [k + 1:end])];
        end
        V = xM;


            function [d_raw, d_prime, r2, r, p] = i_ld(site1, site2)
                n = length(site1);

                stAa = unique(site1);
                stBb = unique(site2);

                Aid = 1;
                Bid = 1;
                fA = sum(site1 == stAa(Aid)) / n;
                fB = sum(site2 == stBb(Bid)) / n;

                %D_raw's sign is arbitrary:
                %A common convention is to set A, B to be the
                %common allele and a, b to be the rare allele

                if (fA < 0.5), fA = 1 - fA;
                    Aid = 2;
                end % Aid tells which one is common allele
                if (fB < 0.5), fB = 1 - fB;
                    Bid = 2;
                end % Bid tells which one is common allele

                x = 0;
                for k = 1:n
                    if (site1(k) == stAa(Aid) && site2(k) == stBb(Bid)),
                        x = x + 1;
                    end
                end
                d_raw = x ./ n - fA * fB;

                if (nargout > 1),
                    fa = 1 - fA;
                    fb = 1 - fB;

                    %r2=(d_raw*d_raw)./(prod(probMajor)*prod(1-probMajor));
                    %r2=min([1 r2]);

                    r2 = (d_raw * d_raw) ./ (fA * fa * fB * fb);
                    r = sqrt(r2);

                    if (d_raw < 0),
                        x = min(fA*fB, fa*fb);
                        r = -1 * r;
                    else
                        x = min(fA*fb, fa*fB);
                    end
                    d_prime = abs(d_raw) ./ x;
                end


                if nargout > 4
                    M = zeros(2);
                    for k = 1:n
                        if (site1(k) == stAa(1) && site2(k) == stBb(1)),
                            M(1, 1) = M(1, 1) + 1;
                        elseif (site1(k) == stAa(1) && site2(k) == stBb(2)),
                            M(1, 2) = M(1, 2) + 1;
                        elseif (site1(k) == stAa(2) && site2(k) == stBb(1)),
                            M(2, 1) = M(2, 1) + 1;
                        elseif (site1(k) == stAa(2) && site2(k) == stBb(2)),
                            M(2, 2) = M(2, 2) + 1;
                        end
                    end
                    [p] = fisherextest(M(1, 1), M(1, 2), M(2, 1), M(2, 2));
                end


                    function [d_raw, d_prime, r2] = i_ld01(site1, site2)
                        n = length(site1);

                        stAa = [0; 1]; %unique(site1);
                        stBb = [0; 1]; %unique(site2);

                        Aid = 1;
                        Bid = 1;
                        fA = sum(site1 == stAa(Aid)) / n;
                        fB = sum(site2 == stBb(Bid)) / n;

                        %D_raw's sign is arbitrary:
                        %A common convention is to set A, B to be the
                        %common allele and a, b to be the rare allele

                        if (fA < 0.5), fA = 1 - fA;
                            Aid = 2;
                        end % Aid tells which one is common allele
                        if (fB < 0.5), fB = 1 - fB;
                            Bid = 2;
                        end % Bid tells which one is common allele

                        x = 0;
                        for (k = 1:n),
                            if (site1(k) == stAa(Aid) && site2(k) == stBb(Bid)),
                                x = x + 1;
                            end
                        end
                        d_raw = x ./ n - fA * fB;

                        if (nargout > 1),
                            fa = 1 - fA;
                            fb = 1 - fB;

                            %r2=(d_raw*d_raw)./(prod(probMajor)*prod(1-probMajor));
                            %r2=min([1 r2]);

                            r2 = (d_raw * d_raw) ./ (fA * fa * fB * fb);
                            r = sqrt(r2);

                            if (d_raw < 0),
                                x = min(fA*fB, fa*fb);
                                r = -1 * r;
                            else
                                x = min(fA*fb, fa*fB);
                            end
                            d_prime = abs(d_raw) ./ x;
                        end
