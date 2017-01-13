function [geno_out] = snp_genotranspose(geno_in,methodid)
%SNP_GENOTRANSPOSE - transposes gneotype matrix

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

%{
geno_in=[1 2 2 2 1 1;...
         3 3 3 4 4 4]
     
gene_out=[1 2 3 3
          2 2 3 4
          1 1 4 4]
%}      

if nargin<2, methodid=2; end

[n,m]=size(geno_in);
geno_out=zeros(m/2,n*2,'uint8');

switch methodid
    case 1
        for j=1:n
            for k=1:2:m
                  x=(k+1)/2;
                  y=j*2-1;
                  geno_out(x,[y y+1]) = geno_in(j,[k k+1]);
            end
        end
    case 2
        for k=1:n
            geno_out(:,[k*2-1 k*2])=reshape([geno_in(k,1:2:end) geno_in(k,2:2:end)],m/2,2);
        end                
end

