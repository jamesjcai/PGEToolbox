function matrixcolor(D,txtlabel,ttxt)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

n=size(D,1);
for i=1:n-1
for j=i+1:n
    D(i,j)=nan;
end
end
for k=1:n
    D(k,k)=nan;
end


    D=[D;zeros(1,size(D,2))];
    D=[D,zeros(size(D,1),1)];
    pcolor(D);
    axis ij
    
axis square
x=copper; x=x(end:-1:1,:);
colormap(x);
colorbar
title(ttxt);

set(gca,'XTick',[1:n]+0.5);
set(gca,'XTickLabel',txtlabel);
set(gca,'YTick',[1:n]+0.5);
set(gca,'YTickLabel',txtlabel);
set(gca,'TickLength',[0 0]);