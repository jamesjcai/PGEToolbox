function [fv]=par_snp_fst(geno1,geno2,methodtag)
%SNP_FST - Fst SNPs from two populations
%
%    [fst]=snp_fst(geno1,geno2,'Weir')
%    [fst]=snp_fst(geno1,geno2,'Wright')
%    [fst]=snp_fst(geno1,geno2,'Hughes')
%
% See also: SNP_FST3

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<3
    methodtag='Weir';
end
if nargin<2
    error('USAGE: [fst]=snp_fst(geno1,geno2,''Weir'');');
end

[m]=snp_marklen(geno1);
fv=zeros(1,m);
genoacell=cell(1,m);
genobcell=cell(1,m);
for k=1:m
    genoacell{k}=geno1(:,2*k-1:2*k);
    genobcell{k}=geno2(:,2*k-1:2*k);
end

parfor k=1:m
    genoa=genoacell{k};
    genob=genobcell{k};
    switch (lower(methodtag))
        case 'weir'
            fv(k)=i_fst_weir_hiploid(genoa,genob);
        case 'wright'
            fv(k)=i_fst_wright(genoa,genob);
        case 'hughes'
            fv(k)=i_fst_hughes05(genoa,genob);
        otherwise
            error('Wrong METHODTAG');
    end    
end

function [f]=i_fst_weir_hiploid(geno1,geno2)
[p1,p2,n1,n2]=i_subpopallelefreq(geno1,geno2);
[f]=fst_weir(n1,n2,p1,p2);
%if p1>p2, f=-f; end
%REF: http://mbe.oxfordjournals.org/cgi/content/full/23/9/1697

function [f]=i_fst_wright(geno1,geno2)
%p1_i = allele freq in pop 1
%p2_i = allele freq in pop 2
%i = number of sites
%REF: http://www.ajhg.org/AJHG/fulltext/S0002-9297(07)63770-7
%but see: http://intl.genome.org/cgi/content/full/15/11/1496 

[p1,p2,~,~,pbar]=i_subpopallelefreq(geno1,geno2);
q1=1-p1; q2=1-p2;
qbar=1-pbar;
if pbar==0 || qbar==0
    f=nan;
else
    f=1-mean([2*p1.*q1, 2*p2.*q2])./(2*pbar.*qbar);
end

function [f]=i_fst_hughes05(geno1,geno2)
%Ref: Hughes et al. 2005  
%p1 and p2 are the frequencies of the first allele in each of two subpopulations
%http://www.genetics.org/cgi/reprint/170/3/1181
[p1,p2]=i_subpopallelefreq(geno1,geno2);
q1=1-p1;
q2=1-p2;
f=1-sum([sqrt(p1.*p2),sqrt(q1.*q2)]);

function [p1,p2,n1,n2,pbar]=i_subpopallelefreq(geno1,geno2)

%s=2;                    % for SNP pairs, num of subpoulations, s=2
n1=sum(geno1(:)~=5);    % hiploid so n1 = n1*2;
n2=sum(geno2(:)~=5);
%n=n1+n2;

%nc = (1/(s-1))*((n1+n2)-(n1^2+n2^2)/(n1+n2));
%nc = (n-(n1^2+n2^2)/n);

[p1,alle1]=snp_maf(geno1);
[p2,alle2]=snp_maf(geno2);
if alle1~=alle2
    p2=1-p2;
end
if nargout>4
[pbar,~]=snp_maf([geno1;geno2]);
end


