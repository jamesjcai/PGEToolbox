function [b] = numvlabel(a)
%NUMVLABEL - Displays integers into vertical format

%in a=[1 2 3 299 3001]
%out b={
%'00003',
%'00020',
%'00090',
%'12391'}


% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


maxlen = length(num2str(max(a)));

n = length(a);
vlen = zeros(1, n);

for k = 1:n
    lz = '';
    astr = num2str(a(k));
    for x = 1:maxlen - length(astr)
        lz = strcat(lz, '0');
    end
    num{k} = strcat(lz, astr);
end

for i = 1:maxlen % each row
    xn = '';
    for j = 1:n % each colonum
        s = num{j};
        xn = strcat(xn, s(i));
    end
    b{i} = xn;
end

if nargout < 1
    for k = 1:maxlen
        disp(b{k})
    end
end
