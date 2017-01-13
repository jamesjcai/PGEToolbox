function [t,tvar] = thetapi(aln,persite,fromsfs)
%THETAPI - Theta from nucleotide diversity
%i.e., Average number of nucleotide differences
%
% Syntax: [t,tvar] = thetapi(aln,jcflag,persite)
%
% t     - Theta_pi estimated from average number of nucleotide differences
% tvar  - Variance of t
% aln   - sequences
% persite - 1 or 0 calculate t per site or per sequence
% fromsfs - 1 or 0 calculate t by using SFS (site-frequency spectrum)

%Calculated here as the sum of 1.0 - sum of site homozygosity accross sites.
%\[ \widehat\theta_\pi=\sum_{i=1}^{i=S}(1-\sum_{j=1}^{j=4}\frac{k_{j,i} \times (k_{j,i}-1)}{n_i \times (n_i - 1)});k_{j,i}>0 \]
%Where $S$ is the number of segregating sites, $k_{j,i}$ is the number of occurences of the $j_{th}$ character state at site $i$, and $n_i$ is the sample size at site $i$. Calculating the statistic in this manner makes it easy to generalize to an arbitrary number of character states per polymorphic site
%Also equivalent to sum of site heterozygosities:
%\[ \widehat\theta_\pi=\frac{\sum_{i=1}^{i=S}2{p_i}{q_i}}{(\genfrac{}{}{0pt}{}{n}{2})} \]
%Also equivalent to mean pairwise differences, but that's slow to calculate.

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if(nargin<3)
    fromsfs=0;
end
if(nargin<2)
    persite=0;
end
if (isstruct(aln)), seq=aln.seq; else seq=aln; end

if (fromsfs)
    [t,tvar] = i_thetapi_sfs(seq);
else
    [t,tvar] = i_thetapi_meannucdiff(seq);
end

if (persite),
    m=size(seq,2);
    t=t/m;
end



function [k,kvar] = i_thetapi_sfs(seq)
%calculate theta_pi from SFS
    [freqall,n] = sfs(seq);
    %freqall looks like [3 2 1 1 1 0 0 0 0];
    nx=1:(n-1);
    k=2*sum((nx.*(n-nx)).*freqall)./(n*(n-1));
    if (nargout>1),
    	kvar=((n+1)*k)/(3*(n-1))+(2*(n^2+n+3)*k^2)/(9*n*(n-1));    % Tajima (1983) eq.
    end



function [k,kvar] = i_thetapi_meannucdiff(seq)
%calculate theta_pi from average number of nucleotide differences
% k (Tajima 1983, equation A3).

    n=size(seq,1);
    % given n DNA seq, we can estimate the scaled mutation rate theta by
    % k = (n!/2!*(n-2)!)^-1 * groupdiffs
    D=0;
    for i=1:n-1
    for j=i+1:n
    	p=sum(seq(i,:)~=seq(j,:));
        D=D+p;
    end
    end
    k=D/(n*(n-1)/2);
    if (nargout>1),
    	kvar=((n+1)*k)/(3*(n-1))+(2*(n^2+n+3)*k^2)/(9*n*(n-1));    % Tajima (1983) eq.
    	%kvar_stochatic=(1/3)*k+(2/9)*k*k;
    	%kvar_sampling = kvar - kvar_stochatic;
    end

