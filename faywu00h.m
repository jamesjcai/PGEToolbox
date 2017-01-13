function [h,hraw] = faywu00h(nsam, Sn, thissfs)
%FAYWU00H - Fay and Wu (2000) H
%  Syntax: [h] = faywu00h(nsam, Sn, thissfs)

%Significantly negative values of H indicate an excess of high-frequency-
%derived alleles, consistent with recent positive selection.

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (Sn==0), h=0; return; end

%Sn=sum(thissfs);
if (length(thissfs)~=nsam-1)
    error('xx');
end

nx=1:(nsam-1);
th=2*sum((nx.^2).*(thissfs))./(nsam*(nsam-1));
tp=2*sum((nx.*(nsam-nx)).*(thissfs))./(nsam*(nsam-1));  %thetapi
tl=sum(nx.*thissfs)./(nsam-1);

%http://www.genetics.org/cgi/content/full/174/3/1431
hraw=tp-th;
%hraw2=tp-tl;

%a1 = sum(1./nx);
%p=thissfs./(Sn*a1);
%sh=sum(2.*p.*p)*(nsam/(nsam-1));

an=sum(1./nx);
bn=sum(1./(nx.^2));
bn2=sum(1./((1:nsam).^2));
t1=Sn./an;

t2=Sn*(Sn-1)./(an.^2+bn);
%t2=t1*t1;   % in Fay's HTEST.C

n=nsam;

%Equation 12
%http://www.genetics.org/cgi/content/full/174/3/1431
%hvar = var(theta_pi - theta_L);
% H=(2*(theta_pi-theta_L))/(2*sqrt(hvar))
hvar=t1*(n-2)/(6*(n-1))+t2*((18*n^2)*(3*n+2)*bn2...
     -(88*n^3+9*n^2-13*n+6))/(9*n*(n-1)^2);
h=hraw./(2*sqrt(hvar));
%h2=hraw2./sqrt(hvar);

