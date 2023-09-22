function [p] = snp_hwpvalue(genodata)
%SNP_HWPVALUE - HW p value
%  Syntax: [p] = snp_hwpvalue(genodata,markinfo)
%
% SNP_HWPVALUE is the Hardy-Weinberg equilibrium p value, which is the
% probability that its deviation from H-W equilibrium could be explained
% by chance.

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

m = size(genodata, 2);
if (mod(m, 2) > 0), error(''); end
locnum = m / 2;
p = zeros(1, locnum);
fid = p;
mid = p;

%fid=markinfo.fatherid;
%mid=markinfo.motherid;
%pedname=markinfo.pedname;
%id=markinfo.id;
%famidx=unique(pedname);
%famnum=length(famidx);

for k = 1:locnum
    locgeno = genodata(:, [k * 2 - 1, k * 2]);
    p(k) = i_hwpvalue(locgeno, fid, mid);
end


    function [p] = i_hwpvalue(locgeno, fatherid, motherid)
        x = i_genomissing(locgeno);
        y = i_notfounderidx(fatherid, motherid);
        locgeno(union(x, y), :) = [];
        het = sum(locgeno(:, 1) ~= locgeno(:, 2));
        nt = unique(locgeno);
        nt(nt == 0) = [];
        homA = sum((locgeno(:, 1) == nt(1)) & (locgeno(:, 1) == locgeno(:, 2)));
        homB = length(locgeno) - (homA + het);
        [p] = hwe(homA, het, homB);

            function [idx] = i_notfounderidx(fid, mid)
                idx = find(((fid == 0) + (mid == 0)) ~= 2); % no Ancestor

                    function [idx] = i_genomissing(locgeno)
                        idx = find(sum(locgeno == 0, 2) > 0);
