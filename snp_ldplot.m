function snp_ldplot(ldinfo,gridon)
%SNP_LDPLOT - LD matrix plot
%
% snp_ldplot(ldinfo)
%
% SEE ALSO: EMLDRUN, SNP_LDPAIR

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

D=ldinfo.dprime'+ldinfo.r2;

if nargin<2
    if size(D,1)<=100
        gridon=1;
    else
        gridon=0;
    end
end


for k=1:size(D,1)
    D(k,k)=nan;
end

if gridon
    D=[D;zeros(1,size(D,2))];
    D=[D,zeros(size(D,1),1)];

    pcolor(D);
    axis ij
else
    imagesc(D);
end
axis square
x=copper; x=x(end:-1:1,:);
colormap(x);
colorbar
title('r^2');
xlabel('D''');


