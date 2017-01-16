function [rx,fu1,fu2,fu,hap1,hap2]=i_ldblock_rallechap(hapthis2,plotit)
    if nargin<1
        plotit=false;
    end

    %ttext=strrep(xfiles(kx,1:end-4),'_','\_');
    [p_maf]=hap_maf(hapthis2);
    
    idx_h=p_maf>=0.1;
    %idx_l=find(p_maf<12/size(hapthis2,1));
    idx_l=find(~idx_h);
    
    hh=hapthis2(:,idx_h);
    [~,sizHap,seqHap,idxHap] = counthaplotype(hh); 

    %x=idxHap==1;
    % sort by freq
         %[~,sortedidx]=sort(idxHap);
         %hapthis3=hapthis2(sortedidx,:);

    % sort by p-distance between distinct haplotypes
        d=dn_pdist(seqHap);
        d1=d(:,1);
        [~,sortedidx]=sort(d1);
        sizHap=sizHap(sortedidx);
        idxHap2=zeros(size(idxHap));
        for kkk=1:length(sortedidx)
            idxHap2(idxHap==sortedidx(kkk))=kkk;
        end      
        [~,sortedidx2]=sort(idxHap2);
        hapthis3=hapthis2(sortedidx2,:);
   
%{        
    figure 
    snp_vhview(hapthis2);
    title(ttext)

    figure 
    snp_vhview(hapthis3);
    hline(sizHap(1)+1)
    hline(sizHap(2)+sizHap(1)+1)
    hline(sizHap(3)+sizHap(1)+sizHap(2)+1)
    title(ttext)
 %} 
    hapthis4=[hapthis3(:,idx_h),hapthis3(:,idx_l)];

    rx=raggedness(hapthis3(:,idx_h));
    fu=fu97fs(hapthis4);
    fu1=fu97fs(hapthis4(1:sizHap(1),:));
    fu2=fu97fs(hapthis4(sizHap(1)+1:end,:));
    
    hap1=hapthis4(1:sizHap(1),:);
    hap2=hapthis4(sizHap(1)+1:end,:);
    
    
if nargout==0||plotit
    %figure 
    ttext='test';
    snp_vhview(hapthis4,false);
    hline(sizHap(1)+1)
    if length(sizHap)>1, hline(sizHap(2)+sizHap(1)+1); end
    if length(sizHap)>2, hline(sizHap(3)+sizHap(1)+sizHap(2)+1); end
    title(ttext)        
    if ~isempty(idx_l)
        %vline(idx_l,'r-')
        %vline(idx_l+1,'r-')
        vline(sum(idx_h)+1,'r-')
    end
    xlabel(sprintf('rag=%f, Fs=%f, Fs1=%f, Fs2=%f',rx,fu,fu1,fu2));    
end


    
