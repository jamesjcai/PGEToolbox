function [Ls,Ln,ssites,nsites] = lsln_hgcorrected(aln,icode)
%LSLN_HGCORRECTED - Estimates Ls and Ln with Hwang-Green (HG) correction
%
%Ls - Number of synonymous sites
%Ln - Number of nonsynonymous sites
%
%REF: Hwang, D. G. & Green, P. (2004) PNAS 101, 13994–14001.

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (nargin<2), icode=1; end
if (isstruct(aln)),
      seq=aln.seq;
else
      seq=aln;
end
%[n,m] = size(seq);

[R]=hgmutrate(seq);
[M]=mutcsqmat(seq);
a=sum(R.*(M==1))./sum(R.*(M>0));
b=sum(R.*(M==2))./sum(R.*(M>0));
a(isnan(a))=0;
b(isnan(b))=0;

ssites=sum(reshape(a,3,length(a)./3));
nsites=sum(reshape(b,3,length(b)./3));

nsites(1)=3-ssites(1);
ssites(end)=3-nsites(end);

%a(isnan(a))=[]; b(isnan(b))=[];
Ls=sum(ssites(:));
Ln=sum(nsites(:));

