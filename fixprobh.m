function [u]=fixprobh(p,h,s,Ne)
%FIXPROBH - fixation probability with dominance factor H
%
%Usage: 
%
% Kimura 1962 - equ (13)
% http://www.genetics.org/cgi/reprint/47/6/713

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin==4, a=2*Ne*s; end
if nargin==3, a=s; end
if a==0, u=p; return; end

    %c=Ne.*s;
    %if h=0.5; fixprobh == fixprob    

    D=2*h-1;
    F = @(x) exp(-2*a.*D.*x.*(1-x)-2*a.*x);
%    F = @(x) exp(-a.*D.*x.*(1-x)-a.*x);
    Q1 = i_quadl(F,p);
    Q2 = quadl(F,0,1);    
    u=Q1./Q2;
    
    
function y=i_quadl(F,x)
    n=length(x);
    y=zeros(1,n);
    for k=1:n 
        y(k) = quadl(F,0,x(k));
    end
    