function [d] = nucdiv2(aln1,aln2,jcflag)
%NUCDIV2 - Nucleotide diversity between two populations
%  Syntax: [d] = nucdiv(aln,jcflag)
%DNAsp : * Nucleotide diversity, Pi (p), the average number of nucleotide differences per site
% between two sequences (Nei 1987, equations 10.5 or 10.6), and its sampling variance (Nei 1987, equation 10.7).
% * Nucleotide diversity, Pi (p), the average number of nucleotide differences per site between
% two sequences (Nei 1987, equations 10.5 or 10.6; see also Nei and Miller 1990), and its
% sampling variance (Nei 1987, equation 10.7). The standard deviation (or standard error) is the square root of the variance.

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $
%
%See also: snp_hapfreqxdist

% In fact: nucdiv = thetapi/total sites

if nargin<3, jcflag=0; end

if (isstruct(aln1)), seq1=aln1.seq; else seq1=aln1; end
if (isstruct(aln2)), seq2=aln2.seq; else seq2=aln2; end

[n1,m1]=size(seq1);
[n2,~]=size(seq2);

	if jcflag
		[d]=i_quickpi2(seq1,seq2,n1,n2,m1,1);
	else
		[d]=i_quickpi2(seq1,seq2,n1,n2,m1,0);
    end






function [d] = i_quickpi2(seq1,seq2,n1,n2,m,jc)
	d=0;
   	ign=0;  % ignored cell where JC distance cannot be computed
	for i=1:n1
	for j=1:n2
	    p=sum(seq1(i,:)~=seq2(j,:))./m;
	    if jc
	      if p>=0.75;
                p=0;
                ign=ign+1;
          else
                p=(-3/4)*log(1-4*p./3);
	      end
        end
	    d=d+p;
	end
    end
	%dv=d/(n*(n-1)/2-ign);
    d=d/(n1*n2-ign);
    
    
    