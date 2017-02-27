function [G]=snp_vhview(haplo,dispfooter)
%SNP_VHVIEW - visual haplotype (VH)
% Syntax: snp_vhview(haplo)
%
%    blue =  common allele
%    yellow = homozygous genotype for the rare allele
%    gray - missing data (N)


% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-04-23 23:03:55 -0500 (Tue, 23 Apr 2013) $
% $LastChangedRevision: 530 $
% $LastChangedBy: jcai $

%Reference:
%Rieder MJ, Taylor SL, Clark AG, Nickerson DA.
%Sequence variation in the human angiotensin converting enzyme.
%Nat Genet. 1999 May;22(1):59-62.
%PMID: 10319862
%
%Website reference:
%http://pga.gs.washington.edu/VH1.html

if nargin<2
    dispfooter=true;
end

G=[];
if isempty(haplo), return; end
if islogical(haplo), haplo=uint8(haplo)+1; end
[n,m]=size(haplo);
G=3*ones(n,m);
for k=1:m
      x=haplo(:,k);
      y=x;
      z=x;
      y(y>4)=[];
      [a,~,c]=unique(y);      
      switch (length(a))
        case (1)
            x(z==a(1))=1;   % Homozygote-Common allele = 1
        case (2)
	      if (sum(c==1)>sum(c==2)),
    		x(z==a(1))=1;   % Homozygote-Common allele = 1
        	x(z==a(2))=2;   % Homozygote-Rare allele = 2
          else
            x(z==a(1))=2;
    		x(z==a(2))=1;
	      end
      end
      G(:,k)=x;
end
G=cat(2,G,ones(n,1)*3);
G=cat(1,G,ones(1,m+1)*3);
G(G>3)=3;

%{
H=[];
for (k=1:2:n)
    H=[H;G(k,:)+G(k+1,:)];
end

H(H==2|H==3)=1;
H(H==4)=3;
%}

pcolor(G);
axis ij;
grid off;

colormap('default');
ax=jet;
colormap(ax([1,40,24],:))

x=axis;
t=max(x(2),x(4));
x(2)=t;
x(4)=t;
axis(x);
axis equal
%axis off
axis tight
xlabel('Markers (SNPs)')
ylabel('Samples (Chromosomes)')
shading flat

if dispfooter
i_dispheader('Visual Haplotype (VH) View')
    disp('Blue   : Homozygote-common allele ')
    disp('Yellow : Homozygote-rare allele')
    disp('Cyan   : Undetermined ')
    disp(' ')
    disp('* Column = Markers (SNPs)')
    disp('** Row = Samples (Chromosomes)')
i_dispfooter
end

