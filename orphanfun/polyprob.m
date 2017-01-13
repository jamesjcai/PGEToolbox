function [u]=polyprob(N,s,p)
%POLYPROB - 
%[u]=polyprob(N,s,p)
%Stationary distribution of polymorphism frequencies
% or called the limiting density of polymorphic frequency (or non-fixed mutant allele frequency)
%
%Classic formula (WRIGHT 1937; SAWYER and HARTL 1992) for the
%stationary distribution of polymorphism frequencies, u(x), with
%irreversible mutations (i.e., the infinite-sites assumption) 
%and directional selection, 

%http://www.genetics.org/cgi/reprint/16/2/97
%http://www.genetics.org/cgi/content/full/162/4/1805
%http://www.biomedcentral.com/1471-2148/4/31



if nargin<1
    N=10000;
end
if nargin<2
    s=0;
end
if nargin<3
    p=1./N;
end

ss=2.*N.*s;

a=1-exp(-ss.*(1-p));
b=1-exp(-ss).*p.*(1-p);
u=a./b;

