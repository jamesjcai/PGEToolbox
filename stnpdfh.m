function [y]=stnpdfh(x,h,s,Ne)
% STNPDFH - stationary distribution of the frequency X of a newly arisen
% mutation under selection
%
% Wright - 1938 - equ ?  (quasi-stationary distribution)
% http://www.pubmedcentral.nih.gov/pagerender.fcgi?artid=1077089
%
% Kimura 1964 - equ (9.20)
% Diffusion models in population genetics 
% http://www.jstor.org/stable/view/3211856?seq=45
%
% Williamson 2004 - equ (1)
% http://www.genetics.org/cgi/reprint/168/1/463
%
% When the mutation is semidominant (i.e., h=0.5), then stnpdfh==stnpdf.
%
% See also: stnpdfsmpl, stnpdf, stnpdfsmpl_test

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin==4, Sx=2*Ne*s; end
if nargin==3, Sx=s; end

    D=1-2*h;
    F = @(u) exp(-4*Sx.*h.*u - 2*Sx.*D.*u.^2);
    
    Q1 = i_quadl(F,x);
    Q2 = quadl(F,0,1);
    
    Q3=exp(4.*Sx.*h.*x+2.*Sx.*D.*x.^2);
    y=2.*Q3./(x.*(1-x)).*(Q1./Q2);
   
   
function y=i_quadl(F,x)
    n=length(x);
    y=zeros(1,n);
    for k=1:n 
        y(k) = quadl(F,x(k),1);
    end

    