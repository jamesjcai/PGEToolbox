function [B, Q] = wallsBQ(S)
%calculate B=B1/(S-1);B1 = the number of pairs of adjacent segregating;S is
%the segregating sites
%sites that are congruent;Q=(B1+|A|)/S A=the set of all distinct partitions
%induced by congruent pairs of segregating sites
%
% author: 'Mary'Xue Yu <yuxueqd@gmail.com>
%

for i = 1:size(S, 2)
    a1 = S(1, i);
    b1 = S(:, i) == a1;
    S(b1, i) = 0;
end
S(S > 0) = 1; % s1=reference sequence,if the same with s1,value=0, if not, value=1;
E = S'; % rows to columns and columns to rows
n = size(E, 1); % number of segregating sites
C = zeros(n-1, 1);
for i = 1:n - 1
    for j = i + 1
        a = E(i, :);
        b = E(j, :);
        C(i, 1) = all(a == b); % compare whether the pairs of adjacent segregating sites
        B1 = sum(C == 1); % are congruent; count the number of congruent sites
    end
end
B = B1 / (n - 1);
[x] = find(C == 1); % get the row number and column umber of congruent sites
Z = E(x, :); % extract the congruent sites rows
A = size(unique(Z, 'rows'), 1); % count the unique congruent sites
Q = (B1 + A) / n;
