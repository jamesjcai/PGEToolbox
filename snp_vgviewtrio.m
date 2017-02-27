function snp_vgviewtrio(geno)
%SNP_VGVIEWTRIO - visual genotype (VG)
% Syntax: snp_vgviewtrio(geno)
%presents complete raw datasets of individuals' genotype data using a display format called a visual genotype (VG)
%The display and interpretation of large genotype data sets can be simplified by using a graphical display. We have found it useful to present complete raw datasets of individuals' genotype data using a display format called a visual genotype (VG) (see Nickerson et al., Nature Genetics, 19:233-240, 1998, and Rieder et al., Nature Genetics, 22:59-60, 1999). This format presents all data in an array of samples (rows) x polymorphic sites (columns) and encodes each diallelic polymorphism according to a general color scheme where:
%
%    blue - homozygous genotype for the common allele
%    red - heterozygous genotype (both common and rare allele)
%    yellow - homozygous genotype for the rare allele
%    gray - missing data (N)
%
%This array format allows one to visually inspect the data across both individual's diplotypes and polymorphic sites to make comparisons.


%For each SNP, blue, yellow, and red boxes indicate whether the individual is homozygous for the
%common allele, heterozygous, or homozygous for the rare allele, respectively.

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-04-23 23:03:55 -0500 (Tue, 23 Apr 2013) $
% $LastChangedRevision: 530 $
% $LastChangedBy: jcai $

%Nickerson DA, Taylor SL, Weiss KM, Clark AG, Hutchinson RG, Stengard J, Salomaa V, Vartiainen E, Boerwinkle E, Sing CF.
%DNA sequence diversity in a 9.7-kb region of the human lipoprotein lipase gene.
%Nat Genet. 1998 Jul;19(3):233-40.
%PMID: 9662394
%
%http://pga.gs.washington.edu/VG2.html

%[idx,geno]=snp_triosort(popid,geno);

G=snp_hhgeno(geno);

%G(find(G==1))=1; % blue = Homozygote-Common allele
%G(find(G==2))=40; % yellow = Homozygote-Rare allele   3
%G(find(G==3))=65; % red = Heterozygote  4
%G(find(G==4))=24; % cyan = missing  2








%G=sortrows(G);
[n,m]=size(G);

Z=nan(n,m*2);


for k=1:2:m*2
Z(1:3:n,k)=G(1:3:n,(k+1)/2);
Z(2:3:n,k)=G(2:3:n,(k+1)/2);
Z(3:3:n,k+1)=G(3:3:n,(k+1)/2);
end



G=Z;
%add extra row and column for pcolor
%don't have to do this for imagesc
G=cat(2,G,ones(n,1)*4);
G=cat(1,G,ones(1,2*m+1)*4);


pcolor(double(G))
colormap('default');
colorindx=jet;
%colormap(colorindx([1,40,64,24],:));    %[40 yellow, 64 red, 24 cyan]
%colormap([0 0 0; 1  1  0;...
%          0.50196 0.50196 0.50196; 0 1 1]);    %[40 yellow, 64 red, 24 cyan]

colormap([0 0 0.5625; 1 1 0; 0.5 0 0; 0.50196 0.50196 0.50196]);
%colormap([blue; yellow; red; gray]);

axis ij
x=axis;
t=max(x(2),x(4));
x(2)=t;
x(4)=t;
axis(x);
axis equal
%axis off
axis tight
xlabel('Markers (SNPs)')
ylabel('Samples (Individuals)')

i_dispheader('Visual Genotype (VG) View')
disp('Blue   : Homozygote-common allele ')
disp('Yellow : Homozygote-rare allele')
disp('Red    : Heterozygote             ')
disp('Cyan   : Undetermined ')
disp(' ')
disp('* Column = Markers (SNPs)')
disp('** Row = Samples (Individuals)')
i_dispfooter
%shading flat
