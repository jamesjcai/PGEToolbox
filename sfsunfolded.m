function [f]=sfsunfolded(seq,ancseq)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (isstruct(seq)), seq=seq.seq; end
[n,m]=size(seq);
[n2,m2]=size(ancseq);
if n2>1
   ancseq=ancseq(1,:); 
   warning('Multiple ANCSEQs given, but we use only the first sequence.');
end
if m~=m2
    error('ANCSEQ is not in the same length with input SEQs.');
end

f=zeros(1,n-1);
for k=1:m
   site=seq(:,k);
   anc=ancseq(k);
   mutn=sum(site~=anc);
   if mutn>0&mutn<n
       f(mutn)=f(mutn)+1;
   end
end


