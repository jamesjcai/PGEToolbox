function [xpehh, coreset] = snp_xpehh(rsid, radius, popcode1, popcode2, showit)
% SNP_XPEHH - XP-EHH scores
% Syntax: snp_xpehh('rs2291725',500000,'JPT+CHB','YRI',true)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


xpehh = nan(1, 2);
if nargin < 5, showit = false; end
if nargin < 4, popcode2 = 'YRI'; end
if nargin < 3, popcode1 = 'JPT+CHB'; end
if (strcmp(popcode1, 'JPT+CHB')), popcode1 = 'JPT%2BCHB'; end
if (strcmp(popcode2, 'JPT+CHB')), popcode2 = 'JPT%2BCHB'; end
if nargin < 2, radius = 1000000; end
if (nargin < 1),
    prompt = {'RSID'};
    def = {'rs2291725'}; % def={'rs996521'};  %def={'rs6060371'}; %def={'rs7662684'};
    dlgTitle = 'Input RSID';
    lineNo = 1;
    answer = inputdlg(prompt, dlgTitle, lineNo, def);
    if ~(isempty(answer)),
        rsid = answer{1};
    else
        return;
    end
end

disp(sprintf('Obtaining chromosomal position of SNP %s...', rsid))
[chrid, pos] = snp_locator(rsid, 1);
if ~(chrid > 0 && pos > 0)
    return;
end

startn = max(1, pos-radius);
endn = min(chrlen(chrid), pos+radius-1);


if chrid > 22
    if chrid == 23
        chrid = 'X';
    elseif chrid == 24
        chrid = 'Y';
    else
        error('Wrong chromosome number.')
    end
    query = sprintf('Chr%s:%d..%d', chrid, startn, endn);
else
    query = sprintf('Chr%d:%d..%d', chrid, startn, endn);
end

fprintf('Dowloading phased haplotype data from %s...\n', query);
fprintf('Length of the region = %d bp\n', endn-startn+1);
[s01, hapldata1, markinfo1] = snp_downloadhaplotype(query, popcode1);
[s02, hapldata2, markinfo2] = snp_downloadhaplotype(query, popcode2);


disp('Extracting common markers....')
[x, a, b] = intersect(markinfo1.pos, markinfo2.pos);
[hapldata1, markinfo1] = i_hap_pickmarker(hapldata1, markinfo1, a);
[hapldata2, markinfo2] = i_hap_pickmarker(hapldata2, markinfo2, b);


disp('Checking for the focal marker....')
    [yes, idx1] = ismember(rsid, markinfo1.rsid);
    if ~yes
        error('Supplied rsid is not found in the region.')
    end
    [yes, idx2] = ismember(rsid, markinfo2.rsid);
    if ~yes
        error('Supplied rsid is not found in the region.')
    end


    xpos = markinfo1.pos;

    if (strcmp(popcode1, 'JPT%2BCHB')), popcode1 = 'JPT+CHB'; end
    if (strcmp(popcode2, 'JPT%2BCHB')), popcode2 = 'JPT+CHB'; end


    fprintf('Calculating EHH for %s in population %s...\n', rsid, popcode1);
        [ehh1, coreset1] = i_snp_ehh(hapldata1, idx1, idx1);
        fprintf('Calculating EHH for %s in population %s...\n', rsid, popcode2);
            [ehh2, coreset2] = i_snp_ehh(hapldata2, idx2, idx2);


            [b1] = i_ehh_leftrightboundries(ehh1);
            [b2] = i_ehh_leftrightboundries(ehh2);

            if (size(ehh1, 1) == 2),
                ids1 = b1(1, 1):b1(1, 2);
                ids2 = b1(2, 1):b1(2, 2);
                ihh1a = trapz(xpos(ids1), ehh1(1, ids1));
                ihh1b = trapz(xpos(ids2), ehh1(2, ids2));
            end
            if (size(ehh2, 1) == 2),
                ids1 = b2(1, 1):b2(1, 2);
                ids2 = b2(2, 1):b2(2, 2);
                ihh2a = trapz(xpos(ids1), ehh2(1, ids1));
                ihh2b = trapz(xpos(ids2), ehh2(2, ids2));
            end

            alleleflipped = false;
            if any(sort(coreset1) ~= sort(coreset2))
                error('xxx');
            else
                if coreset1(1) == coreset2(1)
                    xpehh(1) = log(ihh1a./ihh2a);
                    xpehh(2) = log(ihh1b./ihh2b);
                else
                    alleleflipped = true;
                    xpehh(1) = log(ihh1a./ihh2b);
                    xpehh(2) = log(ihh1b./ihh2a);
                end
            end

            coreset = coreset1;
            X = 'ACGT';


            if (nargout < 1 || showit)

                i_dispheader('XP-EHH')
                fprintf('SNP: %s\n', rsid);
                fprintf('Population: %s vs. %s\n', popcode1, popcode2);
                if any(~isnan(xpehh))
                    if ~alleleflipped
                        fprintf('Allele 1: %s\n', X(coreset(1)));
                        fprintf('Area under haplotype 1 in %s: %f\n', popcode1, ihh1a);
                        fprintf('Area under haplotype 1 in %s: %f\n', popcode2, ihh2a);
                        fprintf('Ratio of two areas: %f\n', abs(ihh1a/ihh2a));
                        fprintf('XP-EHH = %f\n\n', xpehh(1));
                        fprintf('Allele 2: %s\n', X(coreset(2)));
                        fprintf('Area under haplotype 2 in %s: %f\n', popcode1, ihh1b);
                        fprintf('Area under haplotype 2 in %s: %f\n', popcode2, ihh2b);
                        fprintf('Ratio of two areas: %f\n', abs(ihh1b/ihh2b))
                        fprintf('XP-EHH = %f\n\n', xpehh(2));
                    else
                        fprintf('Allele 1: %s\n', X(coreset(1)));
                        fprintf('Area under haplotype 1 in %s: %f\n', popcode1, ihh1a);
                        fprintf('Area under haplotype 1 in %s: %f\n', popcode2, ihh2b);
                        fprintf('Ratio of two areas: %f\n', abs(ihh1a/ihh2b));
                        fprintf('XP-EHH = %f\n\n', xpehh(1));
                        fprintf('Allele 2: %s\n', X(coreset(2)));
                        fprintf('Area under haplotype 2 in %s: %f\n', popcode1, ihh1b);
                        fprintf('Area under haplotype 2 in %s: %f\n', popcode2, ihh2a);
                        fprintf('Ratio of two areas: %f\n', abs(ihh1b/ihh2a))
                        fprintf('XP-EHH = %f\n', xpehh(2));
                    end
                else
                    fprintf('N/A.\n');
                end
                i_dispfooter;

                %xpos=markinfo1.pos;
                %xpos=xpos-min(xpos)+1;
                x1 = ehh1(1, :);
                y1 = ehh1(2, :);
                if ~alleleflipped
                    x2 = ehh2(1, :);
                    y2 = ehh2(2, :);
                else
                    x2 = ehh2(2, :);
                    y2 = ehh2(1, :);
                end
                figure;
                subplot(2, 2, 1)
                plot(xpos, x1, 'r-');
                hold on;
                plot(xpos, x2, 'b-');
                ylabel('EHH')
                xlabel('Chromosomal position (bp)')
                title(sprintf('SNP:%s allele: %s', rsid, X(coreset(1))))
                legend({sprintf('Pop. 1: %s', popcode1), sprintf('Pop. 2: %s', popcode2)})
                hline(0.05, 'g:');

                subplot(2, 2, 2)
                plot(xpos, y1, 'r-');
                hold on;
                plot(xpos, y2, 'b-');
                ylabel('EHH')
                xlabel('Chromosomal position (bp)')
                title(sprintf('SNP:%s allele: %s', rsid, X(coreset(2))))
                legend({sprintf('Pop. 1: %s', popcode1), sprintf('Pop. 2: %s', popcode2)})
                hline(0.05, 'g:');

            end


            function [hap, mrk] = i_hap_pickmarker(hap, mrk, idx)
                if isempty(idx), return; end
                hap = hap(:, idx);
                mrk.rsid = mrk.rsid(idx);
                mrk.pos = mrk.pos(idx);
                mrk.maf = mrk.maf(idx);


                    function [h, coresets] = i_snp_ehh(hapldata, n1, n2)

                        %[maffreq,allemajor,alleminor]=snp_maf(hapldata,1);     % MAF of SNPs

                        [core] = hapldata(:, n1:n2);
                        [coresets, temp, coretype] = unique(core, 'rows');

                        % sort coresets and coretype according to frequencies (descending order)
                        [temp, idx] = sort(grpstats(coretype, coretype, 'length'));
                        idx = idx(end:-1:1);
                        coresets = coresets(idx);
                        newcoretype = zeros(size(coretype));
                        for k = 1:length(idx)
                            newcoretype(coretype == idx(k)) = k;
                        end
                        n = size(coresets, 1);


                        data_left = hapldata(:, 1:n1-1);
                        data_right = hapldata(:, n2+1:end);

                        [temp, m] = size(hapldata);
                        m2 = m - (n2 - n1); % n2-n1 is collasped into core
                        h = ones(n, m2); % initialize result H. In the case of K = 1, let EHH = 1.
                        core_p = zeros(n, 1);

                        for k = 1:n,
                            idx = (newcoretype == k);
                            core_p(k) = mean(idx);
                            dl = data_left(idx, :);
                            dr = data_right(idx, :);
                            for i = 1:n1 - 1,
                                [temp1, y] = counthaplotype(dl(:, i:end));
                                z = sum(y);
                                if z > 1
                                    p = y ./ z;
                                    a = sum(p.^2); % p=relative haplotype frequency
                                    h(k, i) = (a - 1 / z) / (1 - 1 / z);
                                end
                            end
                            for j = n2 + 1:m2,
                                [temp1, y] = counthaplotype(dr(:, 1:j-n1));
                                z = sum(y);
                                if z > 1
                                    p = y ./ z;
                                    a = sum(p.^2);
                                    h(k, j) = (a - 1 / z) / (1 - 1 / z);
                                end
                            end
                        end

                            function [b] = i_ehh_leftrightboundries(ehh)
                                [n, m] = size(ehh);
                                m2 = round(m/2);
                                b = zeros(n, 2);
                                usewarning = false;
                                for k = 1:n
                                    hl = ehh(k, 1:m2);
                                    hr = ehh(k, m2+1:m);
                                    ihl = find(hl <= 0.05);
                                    ihr = find(hr <= 0.05);
                                    if isempty(ihl), b(k, 1) = 1;
                                        usewarning = true;
                                    else b(k, 1) = ihl(end);
                                    end
                                    if isempty(ihr), b(k, 2) = m;
                                        usewarning = true;
                                    else b(k, 2) = m2 + ihr(1);
                                    end
                                end
                                if usewarning
                                    warning('EHH for at least one haplotype does not reach 0.05. Consider longer genomic region.')
                                        %http://www.nature.com/nature/journal/v449/n7164/extref/nature06250-s1.pdf
                                        %page 2: "If, however, EHH doesn't drop in both directions below 0.05 within 2.5MB
                                        %of the core SNP, we skip the iHS test for that SNP."
                                    end
