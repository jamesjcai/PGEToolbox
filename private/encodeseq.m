function [s] = encodeseq(s0,seqtype)
%ENCODESEQ - Converts nucleotide from letters to integer.
%
% Syntax: [s] = encodeseq(s0,seqtype)
%
% Inputs:
%    s0     - Sequences letter representation
%
% Outputs:
%    s    - Sequences integer representation
%
%
% See also: CODONISESEQ, ENCODEALN

% Molecular Biology & Evolution Toolbox, (C) 2007
% Author: James J. Cai
% Email: jamescai@stanford.edu
% Website: http://bioinformatics.org/mbetoolbox/
% Last revision: 5/18/2007
if (nargin<2)
    seqtype=1;
end
s=s0;
switch (seqtype)
    case (1)
         s = i_encode_n(s0);
    case (2)
         s = i_encode_n(s0);
    case (3)
         s = i_encode_a(s0);
end