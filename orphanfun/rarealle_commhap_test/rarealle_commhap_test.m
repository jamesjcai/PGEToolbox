% using prepared haps
popidv = {'yri', 'ceu', 'asn'};

fprintf('[......................]\n');

X2v = [];
c2v = [];
p2v = [];
for kpopid = 2:2
    fprintf('[');
    popid = popidv{kpopid};
    fid = fopen(sprintf('runninglog_%s.txt', popid), 'w');

    for chrid = 18:18
        fprintf('.');
        load(sprintf('../1000genomes/2010_07/genotype/%s/%s_chr%d_geno.mat', ...
            popid, lower(popid), chrid), 'geno');
        load(sprintf('../biodata/1000genomes/2010_07/genotype/%s/%s_chr%d_mark.mat', ...
            popid, lower(popid), chrid), 'mark');
        load(sprintf('%s_chr%d_GABRIELblocks.mat', ...
            popid, chrid), 'b');
        BLK = b;
        clear b;
        for k = 1:min([300, length(BLK)])
            %for k=1:length(BLK)
            lenx = BLK{k}.markers(end) - BLK{k}.markers(1) + 1;
            numx = length(BLK{k}.markers);
            if lenx > 10000 && numx > 10
                [ispos, idx] = ismember(BLK{k}.markers, mark.pos);
                if any(~ispos), error('x'); end
                matfilex = sprintf('hap/%s_chr%d_%d_%d', ...
                    popid, chrid, mark.pos(min(idx)), mark.pos(max(idx)));

                if ~(exist([matfilex, '.mat'], 'file'))
                    picked_i = min(idx):max(idx);
                    [genothis, markthis] = snp_pickmarker(geno, mark, picked_i);
                    [~, ~, ~, isdiallelic] = snp_maf(genothis);
                    [genothis, markthis] = snp_pickmarker(genothis, markthis, isdiallelic);
                    [hapthis2] = i_snp_plemrun(genothis, k, false);
                    save(matfilex, 'hapthis2', 'genothis', 'markthis');
                else
                    load(matfilex, 'hapthis2', 'genothis', 'markthis');
                end
                if size(hapthis2, 1) ~= size(genothis) * 2, continue; end
                [p_maf, ~, minalle] = snp_maf(genothis);

                idx_h = p_maf >= 0.1;
                idx_l = find(p_maf < 2/size(hapthis2, 1));

                if length(idx_l) >= 10
                    minalle_l = minalle(idx_l);

                    hh = hapthis2(:, idx_h);
                    [numHap, sizHap, seqHap, idxHap] = counthaplotype(hh);
                    x = idxHap == 1;
                    a = 0;
                    a2 = 0;
                    c = sizHap(1);
                    d = sum(sizHap) - c;
                    for kk = 1:length(idx_l)
                        y = double(hapthis2(:, idx_l(kk))) == minalle_l(kk);
                        %if ~any(y), error('xx'); end
                        rareidx = find(y);
                        modv = mod(rareidx, 2);

                        %{
                        % swip alleles
                        for kkk=1:length(rareidx)
                            if rand<0.5
                                if modv(kkk)==1
                                    y(rareidx(kkk))=false;
                                    y(rareidx(kkk)+1)=true;
                                else
                                    y(rareidx(kkk))=false;
                                    y(rareidx(kkk)-1)=true;
                                end
                            end
                        end
                        a=a+sum(x&y);
                        a2=a2+sum(~x&y);
                        %}

                        for kkk = 1:length(rareidx)
                            if modv(kkk) == 1
                                isok = (x(rareidx(kkk)+1) && x(rareidx(kkk))) || (~x(rareidx(kkk)+1) && ~x(rareidx(kkk)));
                            else
                                isok = (x(rareidx(kkk)-1) && x(rareidx(kkk))) || (~x(rareidx(kkk)-1) && ~x(rareidx(kkk)));
                            end
                            if isok
                                a = a + sum(x & y);
                                a2 = a2 + sum(~x & y);
                            end
                        end
                    end
                    if a + a2 > 0
                        %[a,a2,c,d]
                        [Px2, ~, X2] = chi2test(a, a2, c, d);
                        if a / a2 < c / d, X2 = -X2; end
                        fprintf(fid, '%d\t%d\t%d\t%d\t%f\t%f\n', ...
                            a, a2, c, d, X2, Px2);
                        X2v = [X2v; X2];
                        c2v = [c2v; c];
                        p2v = [p2v; Px2];
                    end
                end
            end
        end
    end
    fprintf(']\n');
    fclose(fid);
    save(sprintf('%s_test_res', popid), 'X2v', 'c2v', 'p2v');
end
