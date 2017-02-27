function [res,data,colnames] = line4sewfww(aln)
%LINE4SEWFWW - writes input line for SEW or FWW
%dDn=data(:,1);
%dPn=data(:,3);
%dDs=data(:,5);
%dPs=data(:,7);

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


[Ds,Ps,Dn,Pn,Ls,Ln] = dspsdnpn(aln);
[n,m]=size(aln.seq);

haploidratio=1;
colnames = {'Dn' 'Ln(Div)' 'Pn' 'Ln(Poly)' 'Ds' 'Ls(Div)' 'Ps' 'Ls(Poly)' 'n' 'C'};
data=[Dn,Ln,Pn,Ln,Ds,Ls,Ps,Ls,n,1];

res=sprintf('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d',...
Dn,Ln,Pn,Ln,Ds,Ls,Ps,Ls,n,1);


if (nargout<1)
fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',colnames{1:length(colnames)})
disp(res)
end
