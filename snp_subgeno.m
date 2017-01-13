function [geno2]=snp_subgeno(idx,geno,dim)
%[geno2]=snp_subgeno(idx,geno)
%
%See also: SNP_PICKMARKER

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (nargin<3)
	dim=1;        % 1- choose markers; 2 - pick individuals
end

if (dim==1)
    s=idx;
    s2=s*2-1;
    sx=zeros(1,length(s)*2);
    for (k=1:length(s));
        sx(k*2-1)=s2(k);
        sx(k*2)=s2(k)+1;
    end
    geno2 = geno(:,sx);
else
    geno2 = geno(idx,:);
end
%    mark2.rsid=mark.rsid(s);
%    mark2.allele=mark.allele(s);
%    mark2.strand=mark.strand(s);
%    mark2.chr=mark.chr(s);
%    mark2.pos=mark.pos(s);