function [p] = snp_obshet(geno)
%SNP_OBSHET - observed percentage of heterozygosity indiv

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-05-20 19:12:49 -0500 (Mon, 20 May 2013) $
% $LastChangedRevision: 543 $
% $LastChangedBy: jcai $

[n,m]=size(geno);
if(mod(m,2)>0), error(''); end

locnum=m/2;
p=zeros(1,locnum);

for k=1:locnum
      locgeno=geno(:,[k*2-1,k*2]);
      p(k)=i_obshet(locgeno);
end



%~*~*~*~*~*~*~*~*~*
%  SUBFUNCTIONS   *
%~*~*~*~*~*~*~*~*~*

function [p] = i_obshet(locgeno)
locgeno((locgeno(:,1)>4 | locgeno(:,2)>4),:)=[];
% Thank Eran Elhaik for pointing out the bug at version 1.37
locgeno((locgeno(:,1)<1 | locgeno(:,2)<1),:)=[];

%sum(locgeno(:,1)~=locgeno(:,2))
%size(locgeno,1)
p=sum(locgeno(:,1)~=locgeno(:,2))./size(locgeno,1);


function [p] = i_obshet_old(locgeno)
	%[n,m]=size(locgeno); %assert(m==2)
	[numHap,sizHap,seqHap] = counthaplotype(locgeno);
	[n,m]=size(seqHap);
	%assert(m==2);
	het=0;
	hom=0;
	missing=0;
	for (k=1:n),
		allele1=seqHap(k,1);
		allele2=seqHap(k,2);
	      if (allele1~=0 && allele2~=0),
		      if (allele1==allele2),
			hom=hom+sizHap(k);
		      else
			het=het+sizHap(k);
		      end
	      else
		      missing=missing+sizHap(k);
	      end
	end
	%x=sum(sum(geno==0,2)>0)
	p=het/(het+hom);

	% preHet % ref: CheckData.java
	%[a,b]=freqtable(geno);
	%b(find(a==0))=[];
	%a(find(a==0))=[];
	%1-sum(b.^2)/(sum(b)^2)
