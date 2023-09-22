function [t] = time1st2prob(p, p0, N)
%TIME1ST2PROB - Time to reache frequency p for the first time
%
% [t]=time1st2prob(p,p0,N)
% i.e. the mean "first arrival time"
%
% REF: Kimura Ohta 1973 (equ 15, 16 and Appendix)
% http://www.pubmedcentral.nih.gov/picrender.fcgi?artid=1212997&blobtype=pdf
%
% See also: TIME2PROB

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin < 3, N = 10000; end
if nargin < 2, p0 = 0; end

if p0 == 0
    t = 4 .* N .* ((1 - p) .* log(1-p) ./ p + 1);
    %fprintf('4nex=%f\n',2.*N.*p);
else
    t = 4 .* N .* ((1 - p) .* log(1-p) ./ p - (1 - p0) .* log(1-p0) ./ p0);
end
