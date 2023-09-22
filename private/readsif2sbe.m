function [sbeG, sbeNode] = readsif2sbe(filename)

if nargin < 1
    [filename, pathname] = uigetfile( ...
        {'*.txt;*.tab', 'Text Files (*.txt, *.tab)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Pick a Text file');
    if ~(filename), aln = [];
        return;
    end
    filename = [pathname, filename];
end
[a, pp, b] = textread(filename, '%s%s%s');
c = [a; b];

[sbeNode, z, x] = unique(c);
n = numel(sbeNode);
y = reshape(x, [numel(x) / 2, 2]);
%z=[x(1:numel(x)/2),x(numel(x)/2+1:end)];
%all(y==z);

sbeG = false(n);
for k = 1:numel(x) / 2
    sbeG(y(k, 1), y(k, 2)) = true;
    sbeG(y(k, 2), y(k, 1)) = true;
end
