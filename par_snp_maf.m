function [p,allemajor,alleminor] = par_snp_maf(geno,ishaplotype)
%SNP_MAF - Mininum Allele Frequence of SNP
%[p,allemajor,alleminor] = snp_maf(geno)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<2, ishaplotype=false; end
if isempty(geno)
    p=[];allemajor=[];alleminor=[];
    return;
end

if ~ishaplotype
    m2=size(geno,2);
    if(mod(m2,2)>0), error('Wrong GENODATA!'); end
    m=m2/2;
    p=zeros(1,m);
    allemajor=zeros(1,m);
    alleminor=zeros(1,m);
    genocell=cell(1,m);
    for k=1:m
        genocell{k}=geno(:,[k*2-1 k*2]);
    end
    parfor k=1:m
          x=genocell{k};
          x(sum((x==5),2)>0,:)=[]; % remove       

          if isempty(x)
            p(k)=nan;
            allemajor(k)=nan;
            alleminor(k)=nan;
          else
              [a,~,c]=unique(x(:));
              if length(a)==1
                  allemajor(k)=a(1);
                  alleminor(k)=nan;
                  p(k)=0;
              elseif length(a)==2
                  y=sum(c==1)/length(x(:));
                  [y,idx]=min([1-y,y]);
                  p(k)=y;
                  allemajor(k)=a(idx);
                  if y>0 && y<1
                    alleminor(k)=a(sum(idx==1)+1);
                  else
                    alleminor(k)=a(idx);
                  end
              else
                  % deal with more than two alleles
                  c=x(:);
                  d=[sum(c==1),sum(c==2),sum(c==3),sum(c==4)];
                  [~,x2]=max(d);
                  allemajor(k)=x2;
                  d(x2)=0;          % remove largest allele
                  [~,x2]=max(d);   % second largest allele as minor allele
                  y=sum(c==x2)/length(c);
                  p(k)=y;
                  alleminor(k)=x2;
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
[~,m]=size(haplo);
p=zeros(1,m);
alle=5*ones(1,m);
alle2=5*ones(1,m);
for k=1:m
      x=haplo(:,k);
      x(x==5)=[];
      if ~isempty(x)
          [a,~,c]=unique(x);
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


