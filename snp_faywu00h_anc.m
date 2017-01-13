function [hraw]=snp_faywu00h_anc(geno,p,ancp)
%SNP_FAYWU00H_ANC - Fay & Wu's H statistics for SNP
% Syntax: [h]=snp_faywu00h(geno,p)
%
% p   -  derived allele frequency, see SNP_DAF

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<2
    p=[];     % derived allele frequency, see SNP_DAF
end

% ref: http://genome.cshlp.org/content/15/11/1553.long

[smpln,markn]=size(geno);
smpln=smpln*2;

[theh]=i_thetah_anc(geno,[],0,p,0,ancp);
[thew]=snp_thetaw(geno,[],0);
[thepi]=snp_thetapi(geno,[],0);
[Sn]=snp_segsites(geno);
%[h] = faywu00h(smpln, Sn, theh);

hraw=thepi-theh;
h=hraw;




function [theh]=i_thetah_anc(geno,mark,persite,p,showwarning,ancp)
%SNP_THETAH - Fay's theta_H from SNPs
%
% Syntax: [theh]=snp_thetah(geno,mark,persite,p,showwarning)
%
% p -  derived allele frequency, see SNP_DAF

if nargin<5, showwarning=1; end
if nargin<4, p=[]; end
if nargin<3, persite=0; end
if nargin<2, mark=[]; end

[smpln]=snp_samplen(geno);


if isempty(p)
    [p] = snp_maf(geno);
    if showwarning
        disp ('WARNING: SNP_THETAH uses MAF as derived frequency.')
    end
end

theh=sum(2.*sqrt(ancp).*p.^2)*(smpln/(smpln-1));

if (persite),
   [L]=snp_markbplen(mark);
   theh=theh/L;
end
