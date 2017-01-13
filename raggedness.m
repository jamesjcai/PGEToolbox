function [r,d]=raggedness(hap,showhist)
%RAGGEDNESS - 
%
% [r]=raggedness(hap,showhist)
%
%See also: MISMCH
%
%REF: 

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if nargin<2
    showhist=false;
end
    
d=mismch(hap);
x=unique(d);
y=grpstats(d,d,'numel');
y=y/sum(y);

s1=0:max(x);
t1=zeros(1,numel(s1));

[~,idx]=ismember(x,s1);
t1(idx)=y;

if showhist
    %figure;
    bar(s1,t1)
    ylabel('Number')
    xlabel('Nucleotide difference')
    xlim([-1 max(s1)+1]);
end

t2=[t1(2:end),nan];
r=nansum((t2-t1).^2);

%{
if nargout>1
    r2=var(d)./mean(d);
end
if nargout>2
    a=gamfit(d);
    r3=a(1);
end
%}
