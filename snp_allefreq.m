function [minfreq,majfreq,minalle,majalle,mincount,majcount] = snp_allefreq(geno,mark)
%SNP_ALLEFREQ - Allele frequencies of SNP
% Syntax: [minfreq,majfreq,minalle,majalle,mincount,majcount] = snp_allefreq(geno,mark)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-05-22 16:25:15 -0500 (Wed, 22 May 2013) $
% $LastChangedRevision: 544 $
% $LastChangedBy: jcai $


[m2]=size(geno,2);
if(mod(m2,2)>0) 
    error('Wrong GENODATA!'); 
end
m=m2/2;

if nargin<2
    for k=1:m
        rsid{k}=['Mrk_',num2str(k)];
    end
else   
    rsid=mark.rsid;
end


minfreq=zeros(1,m);
majfreq=zeros(1,m);
mincount=zeros(1,m);
majcount=zeros(1,m);
minalle=ones(1,m)*5;
majalle=ones(1,m)*5;

for k=1:2:m2
      x=geno(:,[k k+1]);
      x(sum((x==5),2)>0,:)=[]; % remove
      [a,~,c]=unique(x(:));      
      y=sum(c==1)/length(x(:));      
      idx=(k+1)/2;
      
      minfreq(idx)=min(y,1-y);
      majfreq(idx)=max(y,1-y);
      
      mincount(idx)=min(sum(c==1),sum(c==2));
      majcount(idx)=max(sum(c==1),sum(c==2));
            
      
      if (length(a)==1)
          majalle(idx)=a(1);
      elseif (length(a)==2)
          if y<0.5
            minalle(idx)=a(1);
            majalle(idx)=a(2);
          else
            minalle(idx)=a(2);
            majalle(idx)=a(1);
          end
      end
end


NT='ACGT-';

if (nargout<1),
i_dispheader('Allele Frequencyies')
    	fprintf('#\tMajor freq #\tMinor freq #\trsid\n');
        for k=1:m
	        fprintf(['%d\t%s %.2f %3d\t%s %.2f %3d\t%s\n'], k, NT(majalle(k)), majfreq(k), majcount(k),...
                NT(minalle(k)), minfreq(k), mincount(k),rsid{k});
        end
i_dispfooter
end
