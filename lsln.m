function [Ls, Ln, ssites, nsites] = lsln(aln, ispoly, icode)
%LSLN - Estimates Ls and Ln
%Ls - Number of synonymous sites
%Ln - Number of nonsynonymous sites

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if (nargin < 3), icode = 1; end
if (nargin < 2), ispoly = 1; end

if (isstruct(aln)),
    seq = aln.seq;
else
    seq = aln;
end
[n, m] = size(seq);

if (ispoly)
    [S, N] = getsynnonsynsites(icode);
    %[ns,na] = getsynnonsyndiff(icode);
    [cseq] = codonise64(seq);

    ssites = S(cseq);
    nsites = N(cseq);
    Ls = sum(ssites(:)) ./ n;
    Ln = sum(nsites(:)) ./ n;
else
    [dS, dN, dN_dS, lnL, value] = dc_gy94(aln, 1, 2);
    Ls = value.S;
    Ln = value.N;
end
