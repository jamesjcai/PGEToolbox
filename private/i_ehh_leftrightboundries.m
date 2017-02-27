function [b]=i_ehh_leftrightboundries(ehh,idx)
[n,m]=size(ehh);
%m2=round(m/2);
m2=idx;
b=zeros(n,2);
usewarning=false;
for k=1:n
       hl=ehh(k,1:m2); hr=ehh(k,m2+1:m);
       ihl=find(hl<=0.05); ihr=find(hr<=0.05);
       if isempty(ihl), b(k,1)=1; usewarning=true; else b(k,1)=ihl(end); end
       if isempty(ihr), b(k,2)=m; usewarning=true; else b(k,2)=m2+ihr(1); end
end
if usewarning
    warning('EHH for at least one haplotype does not reach 0.05. Consider longer genomic region.')
%http://www.nature.com/nature/journal/v449/n7164/extref/nature06250-s1.pdf
%page 2: "If, however, EHH doesn't drop in both directions below 0.05 within 2.5MB
%of the core SNP, we skip the iHS test for that SNP."    
end