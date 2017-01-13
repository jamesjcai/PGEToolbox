function [p,c,t] = snp_genofreq(geno,mark)
%SNP_GENOFREQ - returns genotype frequencies
%  Syntax: [p,c,t] = snp_genofreq(geno,mark)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

[m2]=size(geno,2);
m=m2/2;

if nargin<2
    for k=1:m
        rsid{k}=['Mrk_',num2str(k)];
    end
else   
    rsid=mark.rsid;
end

c=zeros(3,m);
t=zeros(3,m);

%%
ddgeno=snp_ddgeno(geno);
%%
for k=1:m
   [x,~,z]=unique(ddgeno(:,k),'rows');
   i55=find(x==55);
   if ~(isempty(i55))
      z(z==i55)=[];
      x(i55)=[];
   end
   %disp(sprintf('------ %d',k))
   for i=1:length(x)
       c(i,k)=sum(z==i);
       t(i,k)=x(i);
   end   
   % sort by frequency
   [c(:,k),idx]=sort(-c(:,k));
   t(:,k)=t(idx,k);
end
c=-c;

p=zeros(3,m);
for k=1:m
    c(:,k)
    t(:,k)
      p(:,k)= c(:,k)./sum(c(:,k));
end  




if (nargout<1),
NT='ACGT-';
ntmap={};
for i=1:4
for j=1:4
    ntmap{i*10+j+1}=[i j];
end
end
ntmap{1}=[5 5];

i_dispheader('Genotyp Frequencies')
    	fprintf('#\tGenotyp freq count\tGenotyp freq count\tGenotyp freq count\trsid\n');
        for k=1:m
	        fprintf('%d\t%s %f %d\t%s %f %d\t%s %f %d\t%s\n',...
                    k, NT(ntmap{1+t(1,k)}), p(1,k),...
                    c(1,k), NT(ntmap{1+t(2,k)}), p(2,k),...
                    c(2,k), NT(ntmap{1+t(3,k)}), p(3,k),...
                    c(3,k),rsid{k});
        end
i_dispfooter
end




