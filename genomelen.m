function [x]=genomelen(autoonly,species)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<1, autoonly=1; end
if nargin<2, species=1; end

x=0;
switch species
    case 1    % human
        for k=1:22, x=x+chrlen(k); end
        if ~autoonly
            x=x+chrlen(23);
            x=x+chrlen(24);            
        end
    otherwise
        error('Unknow.')
end
