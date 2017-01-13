function [f]=fst_weir(n1,n2,p1,p2)
%FST_WEIR - Weir's Fstatistic 
%
% [f]=fst_weir(n1,n2,p1,p2)
% calculated unbiased estimates of FST as described by Weir and Cockerham
% 1984 (see also Weir 1996)
% #  Weir, B.S. 1996. Population substructure. In Genetic data analysis II, pp. 161-173. Sinauer Associates, Sunderland, MA.
% # Weir, B.S. and Cockerham, C.C. 1984. Estimating F-statistics for the
% analysis of population structure. Evolution 38: 1358-1370
%
% Example:
%
%[p1,p2] = meshgrid(0:0.05:1);
%[f]=fst_weir(10,10,p1,p2)
%surf(p1,p2,f)
%xlabel('p_1')
%ylabel('p_2')
%zlabel('F_{ST}')
%
%ref:
% http://mbe.oxfordjournals.org/cgi/content/full/23/9/1697

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

%{
if any(size(p1)>1) || any(size(p2)>1) || any(size(n1)>1) || any(size(n2)>1)
    error('This funciton for single genotype only')
end

if p1==0 && p2==0
    f=0; return;
end
if n1<2 || n2<2
    f=nan; return;
end
%}

s=2;                 % num of subpoulations, s=2
n=n1+n2;

% NC - variance-corrected average sample size
nc = (1/(s-1))*((n1+n2)-(n1.^2+n2.^2)./(n1+n2));
% 
% Weighted frequency 
% weighted average of PA across subpopulations
p_hat=(n1./n).*p1+(n2./n).*p2;



% MSG - mean square error within populations
% MSP - mean square error between populations

MSP=(1/(s-1))*((n1.*(p1-p_hat).^2 + n2.*(p2-p_hat).^2));

%sum([n1-1, n2-1])

MSG=(1./sum([n1-1, n2-1])).*(n1.*p1.*(1-p1)+n2.*p2.*(1-p2));

%if (MSP+(nc-1).*MSG)>0
    f = (MSP-MSG)./(MSP+(nc-1).*MSG);
%else
%    f=nan;
%end

%{
%Ht, Hs
%if nargin<1
%    disp('Example: Ht=0.5; Hs=0.42;')
%    Ht=0.5; Hs=0.42;
%end
%f=(Ht-Hs)/Ht;
%}
