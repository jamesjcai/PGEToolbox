function [y] = hasgap(aln)
%HASGAP - Check if alignment contains gap
%
% Syntax: [y]=hasgap(aln)
%
% Inputs:
%    aln   - Alignment structure
%
% Outputs:
%    y     - 1 or 0
%
%
% See also: REMOVEGAPS

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

% NT = ['A' 'C' 'G' 'T' '-'];
% AA = ['A' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'K' 'L' 'M' 'N' 'P' 'Q' 'R' 'S' 'T' 'V' 'W' 'Y' '*' '-'];
if (isstruct(aln)),
    seq = aln.seq;
    if (aln.seqtype == 1 | aln.seqtype == 2) % DNA/RNA
        y = max(max(seq)) > 4;
    elseif (aln.seqtype == 3) % PROTEIN
        y = max(max(seq)) > 20;
    else
        error('No seqtype.');
    end
else
    seq = aln;
    y = max(seq(:)) > 4;
end
