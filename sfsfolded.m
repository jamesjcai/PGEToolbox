function [f]=sfsfolded(seq)
%SFSFOLDED - 
% [f]=sfsfolded(seq)
%

%we have no chimpanzee sequences and therefore cannot tell the the ancestral allele from the mutant
%Instead of counting mutants,we will count the rarest (sometimes called the minor) allele at each site. This time, however, I will
%omit the invariant sites (1 and 6), which do not contribute to the spectrum.

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if (isstruct(seq)), seq=seq.seq; end
[n,m]=size(seq);

f=zeros(1,n-1);
for k=1:m
   site=seq(:,k);
   [a]=unique(site);
   an=length(a);
   if (an>1)      % ignore 
       allelec=zeros(1,an);
       for x=1:an
           allelec(x)=sum(site==a(x));
       end
       [ancn,idx]=max(allelec);
       anc=site(idx);
       mutn=sum(site~=anc);
       f(mutn)=f(mutn)+1;       
   end   
end
