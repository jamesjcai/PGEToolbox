function [x]=snp_pca(genodata,indvinfo,show3d)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2014-03-26 09:30:26 -0500 (Wed, 26 Mar 2014) $
% $LastChangedRevision: 758 $
% $LastChangedBy: jcai $

if nargin<3, show3d=false; end
if nargin<2||isempty(indvinfo)
    [n]=size(genodata,1);
    indvinfo=mat2cellstr(1:n);
end


G=snp_hhgeno(genodata);
%colx={'m','g','b','c','k'};
%txt=indvinfo;

[~,x]=princomp(double(G));
%figure;
if show3d
    plot3(x(:,1),x(:,2),x(:,3),'o')  % : selects all rows
    zlabel('THIRD AXIS')
else    
    plot(x(:,1),x(:,2),'o')  % : selects all rows
    %{
    hold on
    txtidx=indvinfo.spopid;
    for k=1:5
       %text(x(txtidx==k,1),x(txtidx==k,2),[sprintf('\\color{%s}',colx{k}),txt{k}])  % : selects all rows
       text(x(txtidx==k,1),x(txtidx==k,2),txt{k},'Color',colx{k})  % : selects all rows
    end
    hold off
    %}
    vline(0,'r-')
    hline(0,'r-')
end
    box on;
    grid on
    xlabel('FIRST AXIS')
    ylabel('SECOND AXIS')
    title('Principal Components Analysis')