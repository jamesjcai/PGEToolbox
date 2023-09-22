% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

%h=[0 0.5 1 1.3];
warning off
x = linspace(0, 1, 100);
theta = 20;

gam = -10;
h = [-0.5, -0.3, 0, 0.5, 1];
F = [];
D = [];
for k = 1:length(h)
    D = [D; stnpdfh(x, h(k), gam)];
    [f, fraw] = stnpdfsmpl(15, h(k), gam);
    F = [F; fraw];
end
D = D ./ 2;
F = F ./ 2;

figure;
subplot(2, 2, 1)
plot(x(2:end-1), D(:, 2:end-1))
legend(mat2cellstr(h));
subplot(2, 2, 3)
bar(theta*F')
legend(mat2cellstr(h));

gam = 10;
h = [0, 0.5, 1, 1.3];
F = [];
D = [];
for k = 1:length(h)
    D = [D; stnpdfh(x, h(k), gam)];
    [f, fraw] = stnpdfsmpl(15, h(k), gam);
    F = [F; fraw];
end
D = D ./ 2;
F = F ./ 2;

subplot(2, 2, 2)
plot(x(2:end-1), D(:, 2:end-1))
legend(mat2cellstr(h), 2);
subplot(2, 2, 4)
bar(theta*F')
legend(mat2cellstr(h), 2);

warning on
