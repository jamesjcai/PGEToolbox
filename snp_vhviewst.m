function [haph,hapl]=snp_vhviewst(hap)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

    %ttext=strrep(xfiles(kx,1:end-4),'_','\_');
    [p_maf]=hap_maf(hap);    
    idx_h=p_maf>=0.1;
    idx_l=find(~idx_h);
    
    hh=hap(:,idx_h);
    [~,sizHap,seqHap,idxHap] = counthaplotype(hh); 
    %x=idxHap==1;
    % sort by freq
         %[~,sortedidx]=sort(idxHap);
         %hap3=hap(sortedidx,:);
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
        hap3=hap(sortedidx2,:);
  
        haph=hap3(:,idx_h);
        hapl=hap3(:,idx_l);
        hap4=[hap3(:,idx_h),hap3(:,idx_l)];

    %rx=raggedness(hap3(:,idx_h));
    %fu=fu97fs(hap4);
    %fu1=fu97fs(hap4(1:sizHap(1),:));    
    %fu2=fu97fs(hap4(sizHap(1)+1:end,:));
    
if nargout==0
    %figure 
    ttext='test';
    snp_vhview(hap4,false);
    hline(sizHap(1)+1)
    if length(sizHap)>1, hline(sizHap(2)+sizHap(1)+1); end
    if length(sizHap)>2, hline(sizHap(3)+sizHap(1)+sizHap(2)+1); end
    title(ttext)        
    if ~isempty(idx_l)
        %vline(idx_l,'r-')
        %vline(idx_l+1,'r-')
        vline(sum(idx_h)+1,'r-')
    end
%    xlabel(sprintf('rag=%f, Fs=%f, Fs1=%f, Fs2=%f',rx,fu,fu1,fu2));    
end


    
