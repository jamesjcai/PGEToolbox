function [G]=snp_vgview(geno,ancalle)
%SNP_VGVIEW - visual genotype (VG)
% Syntax: snp_vgview(geno)
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
% $LastChangedDate: 2013-04-28 14:12:08 -0500 (Sun, 28 Apr 2013) $
% $LastChangedRevision: 533 $
% $LastChangedBy: jcai $

%Nickerson DA, Taylor SL, Weiss KM, Clark AG, Hutchinson RG, Stengard J, Salomaa V, Vartiainen E, Boerwinkle E, Sing CF.
%DNA sequence diversity in a 9.7-kb region of the human lipoprotein lipase gene.
%Nat Genet. 1998 Jul;19(3):233-40.
%PMID: 9662394
%
%http://pga.gs.washington.edu/VG2.html

if isempty(geno), return; end
if nargin>1
    G=snp_hhgeno(geno,ancalle);
else
    G=snp_hhgeno(geno);
end
%G(find(G==1))=1; % blue = Homozygote-Common allele
%G(find(G==2))=40; % yellow = Homozygote-Rare allele   3
%G(find(G==3))=65; % red = Heterozygote  4
%G(find(G==4))=24; % cyan = missing  2


%G=sortrows(G);
[n,m]=size(G);

%add extra row and column for pcolor
%don't have to do this for imagesc

if length(G)>401
    imagesc(G);
else
    G2=cat(2,G,ones(n,1)*4);
    G2=cat(1,G2,ones(1,m+1)*4);
    pcolor(double(G2));
    axis ij
end
%colormap('default'); %colorindx=jet;
%colormap(colorindx([1,40,64,24],:));    %[40 yellow, 64 red, 24 cyan]
colormap([0 0 0.5625; 1 1 0; 0.5 0 0; 0.50196 0.50196 0.50196]);

axis equal
%axis off
%shading flat
axis tight
box off
xlabel('Markers (SNPs)')
ylabel('Samples (Individuals)')
opengl software

i_dispheader('Visual Genotype (VG) View')
disp('Blue   : Homozygote-common allele ')
disp('Yellow : Homozygote-rare allele')
disp('Red    : Heterozygote             ')
disp('Gray   : Undetermined ')
disp(' ')
disp('* Column = Markers (SNPs)')
disp('** Row = Samples (Individuals)')
i_dispfooter

