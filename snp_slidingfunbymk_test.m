function [v,pos,pvalue] = snp_slidingfunbymk(funfcn,genodata,winsize,step,varargin)
%SNP_SLIDINGFUNBYMK - Sliding windows analysis by marker
% USAGE: [v,pos] = snp_slidingfunbymk(funfcn,genodata,winsize,step,varargin)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


funfcn = fcnchk(funfcn,length(varargin));
%if ~isa(funfcn,'function_handle'), ...

[n,m2]=size(genodata);
m=m2/2;

%markpos=1:m;
%markstart=min(markpos);
%m=max(markpos)-markstart+1;

if (nargin<4), winsize=5; end
if (nargin<5), step=1; end

if (winsize>m||step>m), error('GENODATA is too short.'); end
if (winsize<step), error('Please select appropriate WINSIZE and STEP'); end
if (winsize<1||step<1), error('Please select appropriate WINSIZE and STEP'); end


sp=1:step:m;
ep=sp+(winsize-1);
ep(find(ep>m))=m;
%sp=sp+markstart-1;
%ep=ep+markstart-1;

nwin=length(sp);
v=zeros(1,nwin);
pvalue=ones(1,nwin);
pos=zeros(1,nwin);


%h = waitbar(0,'Please wait...');
for (k=1:nwin),
     %idx=find(markpos>=sp(k) & markpos <=ep(k));
     idx=sp(k):ep(k);
     pos(k)=sp(k);
     if isempty(idx)
         v(k)=nan;
     else
         [geno]=snp_subgeno(idx,genodata);
         v(k)=feval(funfcn,geno,varargin{:});
         pvalue(k)=i_bootstrpfocal(v(k),idx,m,funfcn,genodata);
         
         if isnan(v(k))
             disp('funfcn returns NaN.')
         end
     end
%     waitbar(k/nwin)
end
%close(h)

    p=1:length(v);
    idx=find(~isnan(v));
    pp=p(idx);
    vv=v(idx);

    xi=find(isnan(v));
    yi=interp1(pp,vv,p(xi));
    v(xi)=yi;


if (nargout<1),
 	%plot(smooth(v,'lowess'));
    %plot(smooth(v,'lowess'));
    %stem(v);
    %barline(v)
    plot(v)


    resavg=mean(v);
	resstd=std(v);
	%score=sum(v>=resavg+resstd)+sum(v<=resavg-resstd);
	%[n,m]=size(res);
	info = ['Data - window size ', num2str(winsize), ' step ', num2str(step)];
	title(info);
	%axis([1 length(dnav) min(dnav)*1.1 max(dnav)*1.1]);
	xlabel('Site Position'); ylabel('Statistic');
    hline(resavg,'r-');
    hline(resavg+resstd,'g:');
    hline(resavg-resstd,'g:');
%	    hold on
%	    plot([1:m],ones(1,m)*resavg,'r')
%	    plot([1:m],ones(1,m)*(resavg+resstd),'g--')
%	    plot([1:m],ones(1,m)*(resavg-resstd),'g--')
%	legend('statistic', 'mean','std',-1) ;
    hold off

end

function p=i_bootstrpfocal(val,idx,m,funfcn,genodata)

bootv=zeros(1,100);
idxok=setdiff(1:m,idx);
n=length(idx);
for k=1:100
  x=randperm(length(idxok));
  idxb=idxok(x(1:n));
  [geno]=snp_subgeno(idxb,genodata);
  bootv(k)=feval(funfcn,geno);  
end
p=sum(bootv<val)/100;






