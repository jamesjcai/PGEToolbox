function [B,Q] = wall99bq(aln)
%WALL99BQ - Wall's B & Q statistic test
%
% Syntax: [B,Q] = wall99bq(aln)
% REF: Wall, J. (1999) Genetical Research 74, pp 65-79

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2014-03-26 09:30:26 -0500 (Wed, 26 Mar 2014) $
% $LastChangedRevision: 758 $
% $LastChangedBy: jcai $

[aln]=extractsegregatingsites(aln,1);
if (isempty(aln)), B=0; Q=realmax; return; end
if (isstruct(aln)), seq=aln.seq; else seq=aln; end

[~,m]=size(seq);  % S=m as the number of segregating sites in the sample.
if (m<2), B=realmax; Q=realmax; return; end

Bprime=0;
switch (nargout)
case (1)
	for k=1:m-1
	      Bprime=Bprime+i_iscongruent(seq(:,k),seq(:,k+1));
    end
	B=Bprime/(m-1);

case (2)

    A={};                  % The set of all distinct partitions induced
                           % by congruent pairs of segregating sites.
    S=0;    
	for k=1:m-1
        [yes,S,A] = i_iscongruent2(seq(:,k),seq(:,k+1),S,A);
        Bprime=Bprime+yes;
    end
           
	B=Bprime/(m-1);
	Q=(Bprime+S)/m;       % there is a typo in the original paper

otherwise
	[B,Q]=wall99bq(aln);
	i_dispheader('Wall''s B & Q Test')	
    fprintf('   B : %f\n',B);
	fprintf('   Q : %f\n',Q);
	i_dispfooter
end




function [yes] = i_iscongruent(a,b)
%Call a pair of segregating
%sites congruent if the subset of the data consisting of
%the two sites contains only two different haplotypes. If
%there has been no recombination between the two
%segregating sites, they will be congruent if and only if
%both mutations lie on the same branch of the unrooted
%tree. Another way of thinking about this is to consider
%each segregating site as an unordered partition of the
%sample, where the subsets correspond to those
%individuals that have the same allele at that particular
%site. (A partition of the sample consists of two disjoint
%subsets whose union is the set of individuals in the
%sample.) Two segregating sites are congruent if and
%only if their corresponding partitions are identical.
    [~,~,z1]=unique(a);
	[~,~,z2]=unique(b);
    z=[z1,z2];
    [numHap] = counthaplotype(z);
    if (numHap==2)
        yes=1;
    else
        yes=0;
    end

function [yes,S,A] = i_iscongruent2(a,b,S,A)
    [~,~,z1]=unique(a);
    [~,~,z2]=unique(b);
    z=[z1,z2];
    [numHap]=counthaplotype(z);
    if (numHap==2)
        yes=1;
        z1a=num2str(z1'); z1b=num2str(3-z1');
        z2a=num2str(z2'); z2b=num2str(3-z2');
        if ~ismember(z1a,A)
                A{length(A)+1}=z1a;
                if ~ismember(z1b,A)
                    A{length(A)+1}=z1a;
                    S=S+1;        % S only count once for a pair of cogruent sites
                end
        end
        if ~ismember(z2a,A), A{length(A)+1}=z2a; end
        if ~ismember(z2b,A), A{length(A)+1}=z2b; end        
    else
        yes=0;
    end
    
    
