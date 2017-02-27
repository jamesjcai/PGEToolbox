function [xgenename,output_var]=genometrack(chrid,regionl,regionr,plotit,loaded_var)
% genometrack(chrid,regionl,regionr)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<4, plotit=false; end
%scalf=1000000;
scalf=1;
if nargin==0
    chrid=1;
    regionl=90000000;
    regionr=95000000;
    plotit=true;
end
if nargin<5
    %load namedgenes_hs_ncbi36_hg18_v7;
    load namedgenes_hs_grch37_hg19_v7
%disp('assembly: Mar. 2006 (NCBI36/hg18)');
%    load namedgenes_hs_grch37_hg19_v7;
%    disp('assembly: Feb. 2009 (GRCh37/hg19)');
else
    xls_chrid=loaded_var{1};
    xls_strand=loaded_var{2};
    xls_genestart=loaded_var{3};
    xls_geneend=loaded_var{4};
    xls_ensembid=loaded_var{5};
    xls_genename=loaded_var{6};
end
if nargout>1
    output_var{1}=xls_chrid;
    output_var{2}=xls_strand;
    output_var{3}=xls_genestart;
    output_var{4}=xls_geneend;
    output_var{5}=xls_ensembid;
    output_var{6}=xls_genename;
end

%xls_chrid	xls_strand	xls_genestart	xls_geneend	xls_ensembid
%xls_genename
idx=(xls_chrid==chrid&xls_genestart>=regionl&xls_geneend<=regionr);
idx2=(xls_chrid==chrid&xls_genestart<regionl&xls_geneend>=regionl);
idx3=(xls_chrid==chrid&xls_genestart<regionr&xls_geneend>=regionr);
idx=idx|idx2|idx3;


xgenename=xls_genename(idx);
xstrand=xls_strand(idx);
xgenestart=xls_genestart(idx);
xgeneend=xls_geneend(idx);
if plotit
if ~isempty(xgenestart)
    i_geneplot([xgenestart xgeneend],xgenename,xstrand,[regionl,regionr],scalf)
else
    plot([regionl,regionr],[0 0],'-k');
end
xlim([regionl,regionr]./scalf)
end

function i_geneplot(params,genex,xstrand,regions,scalf)
    plot(regions./scalf,[0 0],'-k');    
    params=params./scalf;
    nb = size(params,1);  %  number of boxplots
    hw = 0.15;   %  halfwidth of boxes
    %plot([min(min(params)) max(max(params))],[0 0],'k-')
    hold on
    for ii = 1:nb
       temp1 = params(ii,1);
       temp2 = params(ii,2);
       xx = [temp1 temp1 temp2 temp2 temp1];
       % yy = [ii-hw ii+hw ii+hw ii-hw ii-hw];
       yy = [-hw hw hw -hw -hw];
       % plot(xx,yy,'-')
       % fill(xx,yy,'r');
       patch(xx,yy,'w');
       if xstrand(ii)
        plot(temp2,0,'k>');
       else
        plot(temp1,0,'k<');
       end      
       %text(temp1,hw,gene{ii},'fontsize',7,'Rotation',90);
       text(double(temp1+round(0.5*(temp2-temp1))),-hw-(hw*0.3),...
            genex{ii},'fontsize',7,'Rotation',-90);
    end
    hold off
    % make some extra space
    % axlim = axis;
    % axlim(1) = axlim(1)-1;
    % axlim(2) = axlim(2)+1;
    % axlim(3) = axlim(3)-2;
    % axlim(4) = axlim(4)+2;
    % axis(axlim)    
    ylim([-2 1])
    %{
    [a]=get(gca,'position');
    a(2)=0.5;
    a(4)=0.3;
    set(gca,'position',a);
    %}
    set(gca,'yticklabel',[],'ytick',[]);
    set(gca,'xticklabel',[],'xtick',[]);
    %set(gca,'xcolor',get(clf,'color'),'xtick',[]);

