function [d] = nucdiv(aln,jcflag)
%NUCDIV - Nucleotide diversity in a population
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

if (nargin<2),jcflag=0;end

if (isstruct(aln)), seq=aln.seq; else seq=aln; end
[n,m]=size(seq);


if nargout<1
	%Dp=dn_pdist(aln);
	%d1=sum(Dp(:))/(n*(n-1));
	[d1,d2]=i_quickpi(seq,n,m,0);
	reppi=i_bootpi(seq);
	dvar = var(reppi);
        dstd = std(reppi);

	djc=i_quickpi(seq,n,m,1);

	i_dispheader('Nucleotide Diversity (i.e., Nei''s Pi)')
	fprintf('Average number of nucleotide differences: %f\n', d2);
	fprintf('Nucleotide diversity, Pi: %f\n', d1);
	fprintf('    Sampling variance of Pi: %f\n', dvar);
	fprintf('       Standard deviation of Pi: %f\n', dstd);
	fprintf('Nucleotide diversity (Jukes and Cantor), Pi(JC): %f\n', djc);
	i_dispfooter
else
	if (jcflag)
		[d]=i_quickpi(seq,n,m,1);
	else
		[d]=i_quickpi(seq,n,m,0);
	end
end




function [reppi] = i_bootpi(seq,runs,tabs)
	%function simP = boot_Pi(m1,runs,tabs)
	% Does bootstrapping of S = mean Number of pairwise differences between chromosomes
	% tabs sets how often "r" index is printed to screen (default=10)
	% m1 is (n*L) data matrix for Pop1, where n=number chromos in sample and
	% L=No. segregating sites.
	fprintf('[         |          ]\n');
	fprintf(' ');
	if nargin<2, runs=1000; end
	if nargin<3, tabs=50; end
	[n,m]=size(seq);
	reppi=zeros(runs,1);
	for r=1:runs
	  if rem(r,tabs)==0, fprintf(['.']), end		%keep tabs on progress
	  %bootstrap the seq matrix
	  randi = ceil(rand(n,1)*n);
	  reppi(r) = i_quickpi(seq(randi,:),n,m,0);
	end
	fprintf('\n');



function [dv,d] = i_quickpi(seq,n,m,jc)
	d=0; dv=0;
	if (n==1), return; end
	  ign=0;  % ignored cell where JC distance cannot be computed
	for i=1:n-1
	for j=i+1:n
	    p=sum(seq(i,:)~=seq(j,:))./m;
	    if jc
	      if p>=0.75
                p=0;
                ign=ign+1;
          else
                p=(-3/4)*log(1-4*p./3);
	      end
	    end
	    d=d+p;
	end
	end
	dv=d/(n*(n-1)/2-ign);