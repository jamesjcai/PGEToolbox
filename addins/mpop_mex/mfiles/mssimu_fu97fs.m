n = 5000;
c0 = zeros(1, n);
c1 = zeros(1, n);
a0 = zeros(1, n);
a1 = zeros(1, n);
for k = 1:n
    [a, b, c] = fu97fs(ms_exe('120 1 -t 6.4'));
    c0(k) = c;
    a0(k) = a;
end
for k = 1:n
    [a, b, c] = fu97fs(ms_exe('120 1 -t 6.4 -G 6.93 -eG 0.2 0.0 -eN 0.3 0.5'));
    c1(k) = c;
    a1(k) = a;
end

figure;
subplot(2, 2, 1)
histcontour(c0)
subplot(2, 2, 2)
histcontour(c1)
subplot(2, 2, 3)
histcontour(a0)
subplot(2, 2, 4)
histcontour(a1)
