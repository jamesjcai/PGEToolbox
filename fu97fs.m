function [Fs,StrobeckS,Sp] = fu97fs(aln)
%FU97FS - Compute Fu's Fs statistic
%  Syntax: [Fs,StrobeckS] = fu97fs(aln)
%
%The Fs test statistic (Fu 1997, equation 1) is based on the haplotype (gene)
% frequency distribution conditional the value of theta (Ewens 1972, equations 19-21).

% This statistic, which is based on the Ewens' sampling distribution
% (Ewens 1972), has low values with the excess of singleton mutations
% caused by the expansion.

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if (isstruct(aln)), seq=aln.seq; else seq=aln; end

Fs=nan; StrobeckS=nan;
theta=thetapi(seq);
if (theta~=0)
    n=size(seq,1);
	k0=counthaplotype(seq);

	% There is an error in original Fu (97) paper
	Sn=prod(theta+(0:n-1));
   
	% S_k is the coefficient of theta_k in S_n (Ewens 1972; Karlin and Mcgregor 1972)
	%
	%S(n,k) is the coefficient of theata^k in theata(n), itahat is, is the abaolute
	%value of a Stirling number of the first kind.

	X=stirling1(n,n); Stir=X(n,:);
	%Stir=mstirling(n-1); Stir=fliplr(Stir);
	
    k=k0:n;
	Sk=abs(Stir(k));
    Sp=sum(Sk.*theta.^k./Sn);
	Fs=log(Sp/(1-Sp));
    
    if (nargout~=1),
        k=1:k0;
        Sk=abs(Stir(k));
        StrobeckS=sum(Sk.*theta.^k./Sn);
    end
end

if nargout<1
i_dispheader('Fu''s Fs Test and Strobeck''s S Test');
	fprintf ('Fu''s Fs statistic : %f\n',Fs);
	fprintf ('Strobeck''s S statistic : %f\n',StrobeckS);
i_dispfooter;
end
   
%Fu's Fs was negative (Fs = ?92.81; P < 0.0001) which occurs when an excess
%of rare haplotypes is present and suggests that either population expansion 
%or genetic hitchhiking has taken place (Fu 1997).
    
    