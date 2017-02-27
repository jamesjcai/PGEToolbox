popidv={'yri','ceu','asn'};

fprintf('[......................]\n');

for kpopid=2:2
    fprintf('['); 
    popid=popidv{kpopid};
  
for chrid=1:22
    fprintf('.');
    load(sprintf('C:/biodata/1000genomes/2010_07/genotype/output3/%s_chr%d_geno.mat',...
         lower(popid),chrid),'geno');
    load(sprintf('C:/biodata/1000genomes/2010_07/genotype/output3/%s_chr%d_mark.mat',...
         lower(popid),chrid),'mark');
    load(sprintf('C:/biodata/LDBlock/1KG_2010_07/ldblocks/%s_chr%d_GABRIELblocks.mat',...
         popid,chrid),'b');
    BLK=b; clear b;
%for k=1:min([5000 length(BLK)])
for k=1:length(BLK)
    %k
    lenx=BLK{k}.markers(end)-BLK{k}.markers(1)+1;
    numx=length(BLK{k}.markers);
    if lenx>10000 && numx>10 
        %&& lenx<200000
        [ispos,idx]=ismember(BLK{k}.markers, mark.pos);
        %if any(~ispos), error('x'); end

        matfilex=sprintf('hap/%s_chr%d_%d_%d',...
	            popid,chrid,mark.pos(min(idx)),mark.pos(max(idx)));
        %if ~(exist([matfilex,'.mat'],'file'))
            picked_i=min(idx):max(idx);
            [genothis,markthis]=snp_pickmarker(geno,mark,picked_i);
            [p_maf,~,~,isdiallelic]=snp_maf(genothis);
            [genothis,markthis]=snp_pickmarker(genothis,markthis,isdiallelic);
            %[p_maf,~,minalle]=snp_maf(genothis);
        if ~(k==23566 && kpopid==2 && chrid==1)
            [hapthis2]=i_snp_plemrun(genothis,k,false);
            save(matfilex,'hapthis2','genothis','markthis');
        end
        %end
    end
end
end
	fprintf(']\n');
end

