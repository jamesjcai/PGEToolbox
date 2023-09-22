function [p] = cioverlap(ci1, ci2)
%CIOVERLAP - returns the proportion overlap (P) between two CIs
%
% [p]=cioverlap(ci1,ci2)


l1 = 0.5 * (ci1(2) - ci1(1));
l2 = 0.5 * (ci2(2) - ci2(1));

x1 = ci1(1) + l1;
x2 = ci2(1) + l2;

w = mean([l1, l2]); % average error bar height
if x1 >= x2
    d = abs(ci1(1)-ci2(2));
else
    d = abs(ci1(2)-ci2(1));
end

p = d ./ w;
