Rfs=zeros(1,100);
Dse=false(1,100);

seedms_mex;
c=0;
while c<100
c
[segseq,segsites,freqx]=simucode_targetfreq;
%[segseq,segsites]=ms_mex(120,30,0);
%[segseq,segsites,freqx]=simucode;

varstr='68 1 -t 1882 -r 3764 50000 -I 2 24 44 160 -g 1 4303.90 -g 2 41228.0 -eg 0.001070 1 0. -eg 0.0001117 2 0. -en 0.000775 2 9.1e-05 -en 0.000785 2 0.01 -ej 0.001600 2 1';
[segseq,segsites]=ms_exe(varstr);

%%
hap=segseq(1:68,:)+1;
if segsites(end)<1
    segsites=segsites*8000;
end
mark.pos=segsites;
[geno]=snp_hap2geno(hap);
%a=freqx(end)
idx=find(mark.pos==4000);
if isempty(idx)
    idx=round(length(segsites)/2);
end
%b=snp_maf(snp_pickmarker(geno,[],idx))
%figure;
%snp_vgview(geno)
LD=snp_ldblock(geno,mark);

%%
if length(LD)==1&&length(LD(1).markers)>10

    x1=LD(1).markers(1);
    x2=LD(1).markers(end);
    idx2=idx-x1+1;

    %mark.pos(x1)
    %mark.pos(x2)

    genox=snp_pickmarker(geno,[],x1:x2);
    %c=snp_maf(snp_pickmarker(genox,[],idx2))

    %a==b
    %b==c

    hapx=snp_plemrun(genox);
    %d=snp_maf(hapx(:,idx2),1);
    %c==d

    %%
    hapy=[hapx(:,idx2),hapx(:,1:idx2-1),hapx(:,idx2+1:end)];
    [rx,fu1,fu2,fu,hap1,hap2]=i_ldblock_rallechap(hapy,false);

    %if ~any(hap1(:,1)~=2)&&~any(hap2(:,1)~=1)
    %    hline(length(hap1(:,1))+1,'g-');
    %end
    xx=fu97fs_ratio(hap1,hap2);
    if ~isnan(xx)
        c=c+1;
        Rfs(c)=fu97fs_ratio(hap1,hap2);
        if ~any(hap1(:,1)~=2)&&~any(hap2(:,1)~=1)
            Dse(c)=true;
        end
    end
end

end


