function [t] = thetaewens(aln)
%THETAEWENS - Estimates theta (4N*mu) using formula 9.26 in Ewens' book

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (isstruct(aln)), seq = aln.seq;
else seq = aln;
end

[n, m] = size(seq);
[numHap, sizHap, seqHap] = counthaplotype(seq);

neps = 0.00001;
xlow = 0.1;
while (i_kval(xlow, n) > numHap)
    xlow = xlow / 10.0;
end
xhigh = 10.0;
while (i_kval(xhigh, n) < numHap)
    xhigh = xhigh * 10.0;
end

while ((xhigh - xlow) > neps),
    xmid = (xhigh + xlow) / 2.0;
    if (i_kval(xmid, n) > numHap),
        xhigh = xmid;
    else
        xlow = xmid;
    end
end
t = xmid;


    function [k] = i_kval(x, n)
        a = [0:n - 1];
        k = sum(x./(a + x));
