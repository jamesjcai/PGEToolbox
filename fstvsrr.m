%% function fstvsrr
%
% FST and recombnation rate

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

m = 0.0001;
c = 1e-7:1e-3:1e-1;
s = 1e-2;
q=0.5;
p=1e-5;

a=p.^(2.*(c./s));
b=min([1 m.^(c./s)]);

c=(1-b).^2.*a;
d=4-a.*(1+b).^a;

x=c./d