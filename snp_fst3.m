function [f]=snp_fst3(genolist,ntype)
%SNP_FST3 - Fst for genotypes of 3 or more populations
%
% genolist={geno1,geno2,geno3};
% [Mfst]=snp_fst3(genolist,'pairwise')
% [Vfst]=snp_fst3(genolist,'cross')
%
% See also: SNP_FST, FST3, FST

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<2
    ntype='cross';
end

switch ntype
    case 'pairwise'

        n=length(genolist);
        f=zeros(n);
        for i=1:n-1
        for j=i+1:n

            %sum()
            g1=genolist{i};
            g2=genolist{j};
            idx=find(snp_genoprct(g1)>0.5 & snp_genoprct(g2)>0.5);
            [g1]=snp_pickmarker(g1,[],idx);
            [g2]=snp_pickmarker(g2,[],idx);

            f(i,j)=nanmean(snp_fst(g1,g2));
            %fprintf('Fst(%d,%d) = %f\n',i,j,f(i,j));
        end
        end
        
    case 'cross'

        s=length(genolist);
        if s<2
            error('FST need more than one genotype data')
        end
        [nvx,pvx]=i_subpopallelefreq(genolist);
        [n,m2]=size(genolist{1});
        m=m2/2;
        f=zeros(1,m);
        for k=1:m
            if any(isnan(pvx(:,k)))
                idx=~isnan(pvx(:,k));
                f(k)=fst3(nvx(idx,k),pvx(idx,k));
            else
                f(k)=fst3(nvx(:,k),pvx(:,k));
            end
        end
end




function [nv,pv]=i_subpopallelefreq(genolist)

s=length(genolist);
nv=[];
pv=[];
[xmaf,xalle]=snp_maf(genolist{1});


for k=1:s
    geno=genolist{k};
    nv=[nv;2*snp_genoprct(geno).*size(geno,1)];
    [pmaf,palle]=snp_maf(geno);
    idx=find(palle~=xalle);
    if ~isempty(idx)
        pmaf(idx)=1-pmaf(idx);
    end
    pv=[pv;pmaf];
end


