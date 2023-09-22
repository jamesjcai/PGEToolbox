seedms_mex;
h = ms_mex(20, 0, 30);
[g, m] = snp_hap2geno(h+1);
[f, am1, am2] = snp_maf(g);


fid = fopen('ori_hap.txt', 'w');
for k = 1:size(h, 1)
    fprintf(fid, '%d', h(k, :));
    fprintf(fid, '\n');
end
x = sum(h == 1) ./ size(h, 1) >= 0.5;
fprintf(fid, '\n');
fprintf(fid, '%d', x);
fclose(fid);


fid = fopen('ori_geno.txt', 'w');
for k = 1:size(g, 1)
    fprintf(fid, '%d', g(k, :));
    fprintf(fid, '\n');
end
fclose(fid);

snp_writelinkage(g, m, 'input3haploview.phs');

[x, y] = snp_plemrun(g)
[haplodata, hmarkinfo] = snp_phaserun(g, m);

fid = fopen('phase_out.txt', 'w');
for k = 1:size(haplodata, 1)
    fprintf(fid, '%d', haplodata(k, :));
    fprintf(fid, '\n');
end
fclose(fid);
