function y=polymorphicpdf(x,s,Ne)
%POLYMORPHICPDF - Polymorphic site frequency probability density function.
%
% Wright 1938
% Sawyer and Hartle 1992
%
% gam = 2Nes

%x=0.2;
%s=0.001;
%Ne=10000;
%y=polymorphicpdf(x,s,Ne)

gam=2.*Ne.*s;


%if gam~=0
    a=1-exp(-2.*gam.*(1-x));
    b=1-exp(-2.*gam);
    c=x.*(1-x);
    y=a./(b.*c);
    
    y(gam==0)=1./x;
%else
%    y=1./x;
%end



