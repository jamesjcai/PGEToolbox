tic
N = 5000;
%[a,b]=ms_mex(N,4,0);
seedms_mex;
[a, b] = ms_exe(sprintf('%d 1 -t 4', N));

%a=logical(a);
c.N = N;
c.m = 0.001;
c.r = 0.001;
c.s = 0.01;
c.h = 0.5;
c.g = 1000;
c.e = ceil(99999*rand);
[x, y] = mpop_mex(a, b, c);

fprintf('\nmpop -N %d -m %f -r %f -S -s %f -h %f -g %d -i msin -o msout -e %d\n\n', ...
    c.N, c.m, c.r, c.s, c.h, c.g, c.e);

idx = find(y == 0.5);
p = sum(x(:, idx)) ./ length(x(:, 1))

toc
%[geno,mark]=mpopout2genomark(x,y);
%snp_writelinkage(geno,mark,'filename')
%system('java -jar Haploview.jar -pedfile filename.ped -info filename.ped.mrk &')


%{
ld=snp_ldblock(geno,mark);
[genox,markx]=snp_pickmarker(geno,mark,ld.markers(1):ld.markers(end));
hap=snp_plemrun(genox);
%}
