function [ze,zeraw] = zeng06e(nsam, Sn, thissfs)
%ZENG06E - Zeng et al (2006) E
%  Syntax: [ze] = zeng06e(nsam, Sn, thissfs)
%
%Significantly negative values of H indicate an excess of high-frequency-
%derived alleles, consistent with recent positive selection.

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if (Sn==0), ze=0; zeraw=0; return; end

%Sn=sum(thissfs);
if (length(thissfs)~=nsam-1)
    error('xx');
end



nx=1:(nsam-1);
a1 = sum(1./nx);
tw = Sn/a1;

%th=2*sum((nx.^2).*(thissfs))./(nsam*(nsam-1));
%tp=2*sum((nx.*(nsam-nx)).*(thissfs))./(nsam*(nsam-1));  %thetapi

tl=sum(nx.*thissfs)./(nsam-1);


%Equation 8
%http://www.genetics.org/cgi/content/full/174/3/1431
% H:  hraw=tp-th;  i.e., hraw=2*(tp-tl);

zeraw = tl - tw;

%a1 = sum(1./nx);
%p=thissfs./(Sn*a1);
%sh=sum(2.*p.*p)*(nsam/(nsam-1));

an=sum(1./nx);
bn=sum(1./(nx.^2));
bn2=sum(1./([1:nsam].^2));
t1=Sn./an;
t2=Sn*(Sn-1)./(an.^2+bn);

n=nsam;

%Equation 14
%http://www.genetics.org/cgi/content/full/174/3/1431

varzeraw = t1*(n/(2*(n-1))-1/an)+t2*(bn/(an*an)+2*((n/(n-1))^2)*bn...
    - 2*(n*bn-n+1)/((n-1)*an) - (3*n+1)/(n-1));

ze=zeraw/sqrt(varzeraw);
