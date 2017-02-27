function t=coaltime(n)
%COALTIME - Coalescent time 

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

t=zeros(1,length(n-1));
for k=n:-1:2
    %r=1/nchoosek(k,2);
    r=2/(k*(k-1));
    %p=exppdf(x,r);
    t(k)=1/r;
end
t=cumsum(t);
