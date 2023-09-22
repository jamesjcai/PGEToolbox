% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

x = linspace(0, 1, 50);
t = [.1, .2, .3, .4, .5];
x0 = [0.4, 0.5, 0.6];
j = 2;

y = zeros(5, length(x));

for i = 1:length(t)
    for k = 1:length(x)
        y(i, k) = randpdf(x(k), t(i), x0(j));
    end
end

figure;
plot(x, y, '-o')
legend({'t=.1', 't=.2', 't=.3', 't=.4', 't=.5'})
%vline(x0(j))