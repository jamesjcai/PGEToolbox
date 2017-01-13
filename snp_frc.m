function [F]=snp_frc(geno,mark)
%SNP_FRC - fraction of inferred recombinant chromosomes

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

G=snp_hhgeno(geno);
n=snp_marklen(geno);

if (size(G,2)~=n)
    error('something wrong')
end

F=zeros(n);
for i=1:n-1
for j=i+1:n    
    [f1,f2]=i_frc(G(:,i),G(:,j));
    F(i,j)=f1;
    F(j,i)=f2;
end
end



function [f1,f2]=i_frc(gt1,gt2)
gt=[gt1,gt2];
%disp('1   : Homozygote-Common allele ')
%disp('2   : Homozygote-Rare allele')
%disp('3   : Heterozygote')
%disp('4   : Undetermined ')
[idx,x]=find(gt==4);
gt(idx,:)=[];       % remove Undetermined idvs
gt=gt(gt(:,1)<3,:); % remove Heterozygote from core SNP
gt=sortrows(gt,1);

gt1=gt(gt(:,1)==1,2);
gt2=gt(gt(:,1)==2,2);

f1=(2*sum(gt1==2)+sum(gt1==3))/(length(gt1)*2);
f2=(2*sum(gt2==1)+sum(gt2==3))/(length(gt2)*2);






