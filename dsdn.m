function [Ds,Dn] = dsdn(aln1,aln2,icode)
%DSDN - Estimates Ds and Dn
%Ds - Number of synonymous substitutions
%Dn - Number of nonsynonymous substitutions

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if (nargin<3), icode=1; end
if (isstruct(aln1)),
      seq1=aln1.seq;
else
      seq1=aln1;
end

if (isstruct(aln2)),
      seq2=aln2.seq;
else
      seq2=aln2;
end
[~,m1] = size(seq1);
[~,m2] = size(seq2);

if (m1~=m2), error('m1~=m2'); end

[cseq1]=codonise64(seq1);
[cseq2]=codonise64(seq2);
[~,cm1] = size(cseq1);
[~,cm2] = size(cseq2);


[ns,na] = getsynnonsyndiff(icode);

Ds=0; Dn=0;
for k=1:cm1
      csite1=cseq1(:,k);
      csite2=cseq2(:,k);
      if (isempty(intersect(csite1,csite2)))
          [~,~,seqHap1] = counthaplotype(csite1);
          [~,~,seqHap2] = counthaplotype(csite2);
          [cds,cdn] = i_codonsynnonsyn2(seqHap1,seqHap2,ns,na);
          %seqHap=[seqHap1;seqHap2];
          [m_num]=i_mutnum2(decodonise64(seqHap1),decodonise64(seqHap2));

       if ((cds+cdn)==0)
          xds=0; xdn=0;
       else
          xds=(cds*m_num)/(cds+cdn);
          xdn=(cdn*m_num)/(cds+cdn);
       end
           Ds=Ds+xds;
           Dn=Dn+xdn;
      end
end


function [cps,cpn] = i_codonsynnonsyn2(csite1,csite2,ns,na)
	[n1] = size(csite1,1);
	[n2] = size(csite2,1);
	cps=0; cpn=0;

	for i=1:n1
	for j=1:n2
		x=csite1(i); y=csite2(j);
		xx=ns(x,y); yy=na(x,y);
		if (xx>0), cps=cps+xx; end
		if (yy>0), cpn=cpn+yy; end
	end
	end

function [m_num] = i_mutnum2(seq1,seq2)
	m_num=0;
	for k=1:3
		thissite1=seq1(:,k);
		thissite2=seq2(:,k);
		if isempty(intersect(thissite1,thissite2))
		      m_num=m_num+1;
		end
	end
