function [y] = randpdf(x, t, x0)

% gfdrand - gene frequency distribution under random genetic drift
%
% Malecot (1944) obtained the approximation for the probability density for
% large t:
%
%           u(p,x,t) = 6*p*(1-p)*exp(-t/(2*Ne))
%

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


%if nargin<4, Vs=0.0483; end
%if nargin<4, Vs=0.383; end

%if nargin<3, x0=0.5; end
%if nargin<2, t=3; end
%if nargin<1, x=0.4; end

%{
A=1./sqrt(2*pi*Vs.*t);
B=(log(x.*(1-x0)./((1-x).*x0))).^2;
B=exp(-Vs./8-B./(2*Vs.*t));
C=sqrt(x0.*(1-x0))./sqrt((x.*(1-x)).^3);
y=A.*B.*C;
%}

% t = #generations/(2*N)


%q=1-p;

%u(x,t,p)
% Kimura (1955) equ 14.5

%Ci=4*p.*q.*(2*j+1)./(j.*(j+1)).*gegenbauerc(j,1.5,(1-2.*p));
%

r = 1 - 2 * x0;
z = 1 - 2 * x;


j = 1;
Ci = 0;

while (j < 20)
    A = ((2 * j + 1) .* (1 - r.^2)) ./ (j .* (j + 1));
    %A=(2*j+1)./(j.*(j+1));
    B = gegenbauerc(j, 1.5, r);
    C = gegenbauerc(j, 1.5, z);
    D = exp(-0.5*j.*(j + 1).*t);
    incr = A .* B .* C .* D;
    j = j + 1;
    Ci = Ci + incr;
end

%{
while (j<300 && abs(incr) > 1e-5)
    A=((2*j+1).*(1-r.^2))./(j.*(j+1));
    B=gegenbauerc(j,1.5,r);
    C=gegenbauerc(j,1.5,z);
    D=exp(-0.5*j.*(j+1).*t);
    incr = A.*B.*C.*D;
    j=j+1;
    Ci=Ci+incr;
end
%}
y = 4 * x0 .* (1 - x0) .* Ci;
%y=Ci;
