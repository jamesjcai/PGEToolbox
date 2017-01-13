function [p]=snp_hwetest(geno)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

    n=size(geno,2)/2;
    p=-1*ones(1,n);
    for k=1:n
        g=geno(:,[k*2-1 k*2]);
        [a]=unique(g(:));
        if length(a)==2
            %p=sum(c==1)/length(c);
            %if p>0.5, p=1-p; end
            %if p>0.05
            hh=snp_hhgeno(g);
            p(k)=hwe(sum(hh==1),sum(hh==3),sum(hh==2));
        end
    end
end


    