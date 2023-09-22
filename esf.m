function p = esf(nj)

% Ewens?sampling formula, the probability distribution of a configuration
% of alleles in a sample of genes under the infinitely-many-alleles model
% of mutation, is proved by a direct combinatorial argument.

% [8] W.J. Ewens, The sampling theory of selectively neutral alleles,
% Theoret. Popul. Biol. 3 (1972), pp. 87?12. Article | PDF (1165 K) | View
% Record in Scopus | Cited By in Scopus (613)

% The Multivuriate Ewens Distribution (MED), called in genetics the Ewens Sampling
% Formula (ESF), describes a specific probability for the partition of the positive integer
% it into parts. It was discovered by Ewens (1972) as providing the probability of the
% partition of a sample of n selectively equivalent genes into a number of different gene
% types (alleles), either exactly in some models of genetic evolution or as a limiting
% distribution (as the population size becomes indefinitely large) in others.

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-16 11:24:50 -0600 (Wed, 16 Jan 2013) $
% $LastChangedRevision: 371 $
% $LastChangedBy: jcai $


%if nargin<2, theta=1; end

n = sum(nj);
k = length(nj);
x = stirling1(n, n);
Stir = abs(x(n, :));

%p=factorial(n)./(Stir(k)*factorial(k)*prod(nj));
p = prod(k+1:n) ./ (Stir(k) * prod(nj));

%{
i=n:-1:k+1;
if length(i)>=length(nj)
    p1=prod(i(length(nj)+1:length(i))).*prod(i(1:length(nj))./nj)./Stir(k);
else
    p1=prod(i./nj(1:length(i)))*prod(nj(length(i)+1:length(nj)))./Stir(k);
end
%}


%Sn=prod(theta+[0:n-1]);           % Sn(theat)     in equ 3.83
%Prob_k=Stir(k).*theta.^k./Sn;     % equ 3.84 Ewens (2003) book
