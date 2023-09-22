N = 1000;
%N=5;
dip = zeros(1, N);
p = zeros(1, N);
nboot = 500;

seedms_mex;
parfor k = 1:N
    %k
    %hap=ms_mex(120,20)+1;
    %hap=ms_exe('120 1 -t 20 -G 6.93 -eG 0.2 0.0 -eN 0.3 0.5')+1;
    %varstr='68 1 -t 1882 -r 3764 50000 -I 2 24 44 160 -g 1 4303.90 -g 2 41228.0 -eg 0.001070 1 0. -eg 0.0001117 2 0. -en 0.000775 2 9.1e-05 -en 0.000785 2 0.01 -ej 0.001600 2 1';
    varstr = 'structure';
    hap = ms_exe(varstr, k) + 1;
    %hap=ms_exe('120 1 -t 20')+1;
    %figure
    %subplot(1,2,1)
    [~, xm] = raggedness(hap, false);
    %[dip(k), p(k)] = hartigansdipsigniftest(xm, nboot);
    [dip(k)] = hartigansdiptest(xm);
    %title(['dip=',num2str(dip(k),3), ', p=',num2str(p(k),3)])
    %subplot(1,2,2)
    %snp_vhviewst(hap)
end
