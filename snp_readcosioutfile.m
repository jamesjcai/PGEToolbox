function [hapldata,markinfo]=snp_readcosioutfile(hapfile,posfile,removefixed)
%SNP_READCOSIOUTFILE
%
% [h,m]=snp_readcosioutfile('out.hap-1','out.pos-1')

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<3
    removefixed=true;
end
if ~exist(hapfile,'file'), error('SNP_READCOSIOUTFILE:NOFILE','Needs haplotype file'); end
if ~exist(posfile,'file'), error('SNP_READCOSIOUTFILE:NOFILE','Needs position file'); end
    
X=importdata(hapfile);
M=importdata(posfile,'\t',1);
pos=uint32(M.data(:,3))';
hapldata=X(:,3:end);
hapldata(hapldata==2)=0;
hapldata=logical(hapldata);
% (COSI Note: "1" is the derived allele, "2" is the ancestral allele.)
f=sum(hapldata)./size(hapldata,1);
if removefixed
    idx=f>0;
    hapldata=hapldata(:,idx);
    markinfo.daf=f(idx);
    markinfo.rsid=num2cellstr(1:size(hapldata,2));
    markinfo.pos=pos(idx);
else
    markinfo.daf=f;
    markinfo.rsid=num2cellstr(1:size(hapldata,2));
    markinfo.pos=pos;    
end


%{
[g,m]=snp_hap2geno(hap);
figure;
subplot(2,1,1)
snp_vgview(g)
f=snp_maf(g);
subplot(2,1,2)
histsfs(f(f>0))
sum(f==0)./length(f);
%}