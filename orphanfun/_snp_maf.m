function [p,allemajor,alleminor] = snp_maf(geno,ishaplotype)
%SNP_MAF - Mininum Allele Frequence of SNP
%[p,allemajor,alleminor] = snp_maf(geno)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2012-12-09 22:08:33 -0600 (Sun, 09 Dec 2012) $
% $LastChangedRevision: 276 $
% $LastChangedBy: jcai $

if nargin<2
    ishaplotype=0;
end
if isempty(geno)
    p=[];allemajor=[];alleminor=[];
    return;
end

if ~ishaplotype
    [n,m2]=size(geno);
    if(mod(m2,2)>0)
        error('Wrong GENODATA!');
    end
    m=m2/2;
    p=zeros(1,m);
    allemajor=zeros(1,m);
    alleminor=zeros(1,m);

    for k=1:2:m2
          x=geno(:,[k k+1]);
          x(sum((x==5),2)>0,:)=[]; % remove
          cidx=(k+1)/2;

          if isempty(x)
            p(cidx)=nan;
            allemajor(cidx)=5;
            alleminor(cidx)=5;
          else
            [a,b,c]=unique(x(:));
              if length(a)<3
                  y=sum(c==1)/length(x(:));
                  [y,idx]=min([1-y,y]);
                  p(cidx)=y;
                  allemajor(cidx)=a(idx);
                  if y>0 && y<1
                    alleminor(cidx)=a(sum(idx==1)+1);
                  else
                    alleminor(cidx)=a(idx);
                  end
              else
                  % deal with more than two alleles
                  c=x(:);
                  d=[sum(c==1),sum(c==2),sum(c==3),sum(c==4)];
                  [x1,x2]=max(d);
                  allemajor(cidx)=x2;

                  d(x2)=0;          % remove largest allele
                  [x1,x2]=max(d);   % second largest allele as minor allele
                  y=sum(c==x2)/length(c);
                  p(cidx)=y;
                  alleminor(cidx)=x2;
              end
          end
    end

else
    haplodata=geno;
    [p,allemajor,alleminor]=i_hapmaf(haplodata);
end



function [p,alle,alle2]=i_hapmaf(haplo)
%HAP_MAF - Mininum Allele Frequence of haplotype
%[p] = hap_maf(haplodata)
%
%n - individuls; m - marker number
[n,m]=size(haplo);
p=zeros(1,m);
alle=5*ones(1,m);
alle2=5*ones(1,m);
for (k=1:m),
      x=haplo(:,k);
      x(x==5)=[];
      if ~isempty(x)
          [a,temp,c]=unique(x);
          y=sum(c==1)/length(x);
          [y,idx]=min([1-y,y]);
          p(k)=y;
          alle(k)=a(idx);
          b=setdiff(a,a(idx));
          alle2(k)=b(1);
      else
          p(k)=nan;
     %     alle(k)=nan;
     %     alle2(k)=nan;
      end
end



%{
function [p,alle,alle2]=i_hapmaf(haplo)
%HAP_MAF - Mininum Allele Frequence of haplotype
%[p] = hap_maf(haplodata)
%
%n - individuls; m - marker number
[n,m]=size(haplo);
p=zeros(1,m);
alle=5*ones(1,m);
for (k=1:m),
      x=haplo(:,k);
      x(find(sum((x==5),2)>0),:)=[]; % remove
      if any(x)
          [a,b,c]=unique(x(:));
          y=sum(c==1)/length(x(:));
          [y,idx]=min([1-y,y]);
          p(k)=y;
          alle(k)=a(idx);
      end
end
alle2=alle;
%}


