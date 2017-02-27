function [p,majalle,minalle,isdiallelic] = snp_maf(geno,ishaplotype)
%SNP_MAF - Mininum Allele Frequence of SNP
%[p,majalle,minalle] = snp_maf(geno)
%
% NOTE: When p=0.5, the major allele will be assigned to first allele
% accoring to the alphabetic order of A, C, G and T. 

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<2, ishaplotype=false; end
if isempty(geno)
    p=[];majalle=[];minalle=[];
    return;
end

if ~ishaplotype
    m2=size(geno,2);
    if(mod(m2,2)>0), error('Wrong GENODATA!'); end
    m=m2/2;
    p=zeros(1,m);
    majalle=zeros(1,m);
    minalle=zeros(1,m);
    isdiallelic=true(1,m);
    for k=1:2:m2
          x=geno(:,[k k+1]);
          x(sum((x==5),2)>0,:)=[]; % remove
          cidx=(k+1)/2;

          if isempty(x)
            p(cidx)=nan;
            majalle(cidx)=nan;
            minalle(cidx)=nan;
            isdiallelic(cidx)=false;
          else
              [a,~,c]=unique(x(:));
              if length(a)==1
                  majalle(cidx)=a(1);
                  minalle(cidx)=nan;
                  p(cidx)=0;
                  isdiallelic(cidx)=false;                  
              elseif length(a)==2
                  y=sum(c==1)/length(x(:));
                  [y,idx]=min([1-y,y]);
                  p(cidx)=y;
                  majalle(cidx)=a(idx);
                  if y>0 && y<1
                    minalle(cidx)=a(sum(idx==1)+1);
                  else
                    minalle(cidx)=a(idx);
                  end
              else
                  % deal with more than two alleles
                  c=x(:);
                  d=[sum(c==1),sum(c==2),sum(c==3),sum(c==4)];
                  [~,x2]=max(d);
                  majalle(cidx)=x2;
                  d(x2)=0;          % remove largest allele
                  [~,x2]=max(d);   % second largest allele as minor allele
                  y=sum(c==x2)/length(c);
                  p(cidx)=y;
                  minalle(cidx)=x2;
                  isdiallelic(cidx)=false;                  
              end
          end
    end

else
    haplodata=geno;
    [p,majalle,minalle]=i_hapmaf(haplodata);
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
          if length(a)==1
              p(k)=0; alle(k)=a; alle2(k)=a;
          else
              y=sum(c==1)/length(x);
              [y,idx]=min([1-y,y]);
              p(k)=y;
              alle(k)=a(idx);
              b=setdiff(a,a(idx));
              alle2(k)=b(1);
          end
      else
          p(k)=nan;
            % alle(k)=nan;
            % alle2(k)=nan;
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


