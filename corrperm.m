function [lowci,upperci,corrsamp,corrobs]=corrperm(x,y,m)
%CORRPERM - permutation test for the correlation between two variables;
%
% [lowci,upperci,corrsamp,corrobs]=corrperm(x,y,m)
%
% We will fix y and permute x;
% Then calculate the correlation between x & y;
% Calcuate the m correlations;

%SOURCE: http://www-stat.stanford.edu/~susan/phylo/index/node63.html

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

n=length(x);
corrobs=corrcoef(x,y);
corrsamp=zeros(1,m);

for i=1:m,
    ind=randperm(n);
    x=x(ind);
    r=corrcoef(x,y);
    corrsamp(i)=r(2,1);
end;

% Get a 95% confidence intervals;

[sortcorr,ignore]=sort(corrsamp);
lowci=sortcorr(0.975*m);
upperci=sortcorr(0.025*m+1);



% xy=[
% 576 3.39;
% 635 3.30;
% 558 2.81;
% 578 3.03; 
% 666 3.44;
% 580 3.07;
% 555 3.00;
% 661 3.43;
% 651 3.36;
% 605 3.13;
% 653 3.12;
% 575 2.74;
% 545 2.76;
% 572 2.88;
% 594 2.96];
% x=xy(:,1);
% y=xy(:,2);
% x'*y
% [lowci,upperci,corrs]=corrperm(x,y,1000)
% hist(corrs,30);
% title('Permutaion test for correlation');
% xlabel('correlation of the x & y');

