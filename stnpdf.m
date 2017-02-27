function [y]=stnpdf(x,s,Ne)
% STNPDF - stationary distribution of the frequency P of a newly arisen
% mutation under selection
% [y]=stnpdf(x,s,Ne)


%[Fisher 1930, Wright, 1969]
%
% Wright 1937 - 
% http://www.pnas.org/cgi/reprint/23/6/307
%
%
% Wright 1938 - equ 39
% http://www.pubmedcentral.nih.gov/pagerender.fcgi?artid=1077089
%
% Kimura 1962 - equ (9.21)
% http://www.jstor.org/stable/view/3211856?seq=45
%
% equ (1)
% http://www.soe.ucsc.edu/research/compbio/papers/ultras_meth_v31.pdf

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin==3, Sx=2*Ne*s; end
if nargin==2, Sx=s; end

Q1=1-exp(-2*Sx.*(1-x));
Q2=1-exp(-2*Sx);
Q3=2./(x.*(1-x));

y=Q3.*(Q1./Q2);

%2*(1-exp(-a.*(1-p)))./((1-exp(-a)).*(p.*(1-p)))


