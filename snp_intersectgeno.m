function [geno1,mark1,geno2,mark2]=snp_intersectgeno(geno1,mark1,geno2,mark2)
%SNP_INTERSECTGENO - returns genotyped data common to both GENO1 and GENO2.
%disp('Extracting common markers....')

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<1
    disp('Reading GENODATA from population 1....')
       [geno1,mark1] = snp_readhapmap;
       
    disp('Reading GENODATA from population 2....')
       [geno2,mark2] = snp_readhapmap; 
end

geno=[]; mark=[];

if ~(isempty(geno1)||isempty(geno2))

x=intersect(mark1.pos,mark2.pos);
[a]=find(ismember(mark1.pos,x));
[b]=find(ismember(mark2.pos,x));
[geno1,mark1]=snp_pickmarker(geno1,mark1,a);
[geno2,mark2]=snp_pickmarker(geno2,mark2,b);

%n1=snp_samplen(geno1)/2;
%n2=snp_samplen(geno2)/2;
%geno=[geno1;geno2];
%mark=mark1;
%mark.population=[ones(n1,1);2*ones(n2,1)];

end
