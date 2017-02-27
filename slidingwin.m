function [res,resavg,resstd,score]=slidingwin(d,winsize,stepsize,showit)
%SLIDINGWIN - Sliding windows anlaysis helper function

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if (nargin<3),
    stepsize=1;
end
if (nargin<2),
    winsize=30;
end

[n,m]=size(d);
res = slidingavg(d,winsize);
resavg=mean(res);
resstd=std(res);
score=sum(res>=resavg+resstd)+sum(res<=resavg-resstd);

if (showit==1),
	plot(res);
	info = ['Data - window size ', num2str(winsize), ' score=', num2str(score)];
	title(info);
	%axis([1 length(dnav) min(dnav)*1.1 max(dnav)*1.1]);
	xlabel('Base (bp)'); ylabel('Data');
    %hold on
    %plot([1:m],ones(1,m)*resavg,'r')
    %plot([1:m],ones(1,m)*(resavg+resstd),'g--')
    %plot([1:m],ones(1,m)*(resavg-resstd),'g--')
	%legend('v', 'mean','std') ;
    %hold off
end

