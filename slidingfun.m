function [v] = slidingfun(funfcn,aln,winsize,stepsize,varargin)
%SLIDINGFUN - Sliding windows analysis
%   USAGE: [v]=slidingfun(funfcn,aln,winsize,stepsize,varargin)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (isstruct(aln)), seq=aln.seq; else seq=aln; end

funfcn = fcnchk(funfcn,length(varargin));
%if ~isa(funfcn,'function_handle'), ...
if (nargin<3), winsize=20; end
if (nargin<4), stepsize=1; end

[n,m]=size(seq);
if (winsize>m||stepsize>m), error('Sequence too short.'); end
if (winsize<stepsize), error('Please select appropriate winsize and stepsize'); end
if (winsize<1||stepsize<1), error('Please select appropriate winsize and stepsize'); end

sp = 1:stepsize:m;
ep = sp + (winsize-1);
ep(ep>m)=m;
nwin = length(sp)-1;
v=zeros(1,nwin);


%h = waitbar(0,'Please wait...');
for k=1:nwin
    seq(:,sp(k):ep(k))
      v(k)=feval(funfcn,seq(:,sp(k):ep(k)),varargin{:});
      % waitbar(k/nwin)
end
% close(h)
if nargout<1
	plot(v);
	res=v;
	resavg=mean(res);
	resstd=std(res);
	scorex=sum(res>=resavg+resstd)+sum(res<=resavg-resstd);
	[n,m]=size(res);
	info = ['Data - window size ', num2str(winsize), ' stepsize ', num2str(stepsize)];
	title(info);
	%axis([1 length(dnav) min(dnav)*1.1 max(dnav)*1.1]);
	xlabel('Site Position'); ylabel('Statistic');
	%    hold on
	%    plot([1:m],ones(1,m)*resavg,'r')
	%    plot([1:m],ones(1,m)*(resavg+resstd),'g--')
	%    plot([1:m],ones(1,m)*(resavg-resstd),'g--')
	%legend('statistic', 'mean','std',-1) ;
    %hold off
end
