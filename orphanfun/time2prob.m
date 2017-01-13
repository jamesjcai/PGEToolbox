function [t]=time2prob(p,N)
%TIME2PROB - Age of a mutant allele with frequency P.
%
%[t]=time2prob(p,N)
%
% To understand the nature of extant variation, we estimated the age of the
% mutant allele within populations in terms of how many generations the
% mutant allele has persisted in the populations since it appeared by
% mutation.
%
% REF: Kimura Ohta 1973 (equ 13)
% http://www.pubmedcentral.nih.gov/picrender.fcgi?artid=1212997&blobtype=pdf
% Applications e.g., Stephens JC (1998) CCR5-delta32
% Book: Mathematical Population Genetics, WJ Ewens, (equ. 3.4 and 3.11)
% 
% See also: TIME1ST2PROB

if nargin<2, N=10000; end
t=-4.*N.*(p.*log(p)./(1-p));             % (equ 13)

