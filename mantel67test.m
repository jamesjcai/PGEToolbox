function [zm,zmv]=mantel67test(a,b)
% mantel67test - Mantel test between two distance matrices
%
% [d]=mantel67test(M,N)
%
%Mantel test (correlation between two distance matrices (in C).)
%REF: Mantel, N. (1967) The detection of disease clustering and a
%generalized regression approach. Cancer Research, 27, 209–220.

% http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/mantel.html

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

[n,m]=size(a);


a=a-tril(a);
b=b-tril(b);
%The original form of the Mantel statistic, zM, is simply the sum of cross
%products (inner product) of unstandardized values in two proximity
%matrices unfolded into two vectors, excluding the trivial diagonal values
%(Mantel, 1967)
zm=sum(sum(a.*b))/(n-1);


%A second, and more general approach is suggested where the proximity
%values in each of the two vectors are standardized before computing the
%Mantel statistic (Mantel, 1967; Mantel & Valand, 1970; Hubert, 1978;
%Sokal, 1979). The standardized Mantel statistic is analogous to r and is
%called rM
%rm=(1/((n-1)-1))*sum(sum(zscore(a).*zscore(b)));


zmv=zeros(1,1000);
for k=1:1000
    nr=randperm(n);
    newa=zeros(n);
    
    for i=1:n-1
    for j=i+1:n
        ni=nr(i);
        nj=nr(j);
        if ni>nj
            temp=nj; nj=ni; ni=temp;
        end
        newa(i,j)=a(ni,nj);
    end
    end
    

%   newa = a(nr,:);          % permute rows
%   newa = newa(:,nr);       % permute cols
   
   
    
zmv(k)=sum(sum(newa.*b))/(n-1);
end

if nargout==0
    hist(zmv,30)
    vline(zm);
end






%http://www.db.dk/binaries/JASIST_%20part2_preprint.pdf
%The Mantel test is sometimes referred to as the quadratic assignment
%procedure (QAP) (Hubert & Schultz, 1976; Hubert, 1987). Quadratic
%assignment is a mathematical problem in permutations.


