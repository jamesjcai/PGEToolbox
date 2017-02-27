function hap_coloredsqu(hap,hapfreq)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

%{
hap=[ 3    3    3    3    4    3    4    1    3    4    3
    3    3    3    3    4    3    4    1    4    4    3
    1    1    1    3    4    1    4    1    3    4    3
    3    3    3    1    4    3    4    1    3    4    3
    1    1    3    3    2    3    2    3    3    2    1
    3    1    1    3    4    1    4    1    3    4    3
    1    3    3    3    2    3    2    3    3    2    1];

hapfreq=[    0.5850
    0.1170
    0.0980
    0.0900
    0.0750
    0.0170
    0.0170];
%}

[n,m]=size(hap);
s=hap(1,:);
G=zeros(size(hap));
for k=1:size(hap,1)
    G(k,:)=hap(k,:)==s;
end

G=cat(2,G,ones(n,1)*3);
G=cat(1,G,ones(1,m+1)*3);
G(G>3)=3;
pcolor(G);
axis ij;
grid off;

colormap('default');
ax=jet;
colormap(ax([64,1,24],:))

%colormap(colorindx([1,40,64,24],:));    %[40 yellow, 64 red, 24 cyan]
%colormap([0 0 0.5625; 1 1 0; 0.5 0 0; 0.50196 0.50196 0.50196]);


x=axis;
t=max(x(2),x(4));
x(2)=t;
x(4)=t;
axis(x);
axis equal
%axis off
axis tight
xlabel('Marker')
ylabel('Haplotype')
for k=1:n
    text(m+1.2,k+0.5,sprintf('%.2f',hapfreq(k)));
end
box off


