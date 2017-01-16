function [Rfs] = fu97fs_ratio(aln1,aln2)
%FU97FS - Compute Fu's Fs statistic
%  Syntax: [Fs,StrobeckS] = fu97fs(aln)
%
%The Fs test statistic (Fu 1997, equation 1) is based on the haplotype (gene)
% frequency distribution conditional the value of theta (Ewens 1972, equations 19-21).

% This statistic, which is based on the Ewens' sampling distribution
% (Ewens 1972), has low values with the excess of singleton mutations
% caused by the expansion.

if (isstruct(aln1)), seq1=aln1.seq; else seq1=aln1; end
if (isstruct(aln2)), seq2=aln2.seq; else seq2=aln2; end

[~,Sp1]=i_fu97fs(seq1);
[~,Sp2]=i_fu97fs(seq2);

Rfs=log(Sp1/Sp2);




function [Fs,Sp]=i_fu97fs(seq)
    Fs=nan; Sp=nan;
theta=thetapi(seq);
if (theta~=0)
    n=size(seq,1);
	k0=counthaplotype(seq);
	Sn=prod(theta+[0:n-1]);
	X=stirling1(n,n); Stir=X(n,:);
    k=k0:n;
	Sk=abs(Stir(k));
    Sp=sum(Sk.*theta.^k./Sn);
	Fs=log(Sp/(1-Sp));
end
    
    