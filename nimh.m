function [ni_mh]=nimh(dn,ds,pn,ps)
%NI by using the Mantel-Haenszel procedure (Adam Eyre-Walker)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if ~(size(dn,2)==1&&size(ds,2)==1&&size(pn,2)==1&&size(ps,2)==1)
    error('xxx')
end

%{
fprintf('NI_{true} = %f\n',sum(ds.*pn)/sum(dn.*ps));
dn=poissrnd(dn,1000,1);
ds=poissrnd(ds,1000,1);
pn=poissrnd(pn,1000,1);
ps=poissrnd(ps,1000,1);
%}

z=dn+ds+pn+ps;
ni_mh=(nansum(ds.*pn./z))/(nansum(dn.*ps./z));

%fprintf('NI_{MH} = %f\n', ni_mh);
%fprintf('NI_{1} = %f\n', sum(ds.*pn)./sum(dn.*ps));
%fprintf('NI_{2} = %f\n', (sum(ds).*sum(pn))./(sum(dn).*sum(ps)));
