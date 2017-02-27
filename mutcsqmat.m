function [MM]=mutcsqmat(seq,icode)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if nargin<2, icode=1; end

[ct,CODON] = codontable(icode);
stops=find(ct=='*');

[n,m]=size(seq);

MM=[];

for ss=1:n
M=zeros(4,m);
for k=1:3:m-2
    %a=seq(k); b=seq(k+1); c=seq(k+2);
    x=seq(ss,k:k+2);
    xc=codonise64(x);
    if ~ismember(xc,stops)
    for i=1:4
        for j=1:3
            y=x;            
            y(j)=i;
            yc=codonise64(y);
            if yc~=xc && xc<65 && yc<65 && ~ismember(yc,stops)
                if ct(yc)==ct(xc)
                    M(i,k+j-1)=1;
                else
                    M(i,k+j-1)=2;
                end
            end
        end        
    end
    end
end
MM=[MM;M];
end