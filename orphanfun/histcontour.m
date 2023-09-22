function histcontour(x, nbin)

if nargin < 2
    nbin = 50;
end
n = length(x);

[count, xout] = hist(x, nbin);
dx = xout(2) - xout(1);
f = count / n / dx; % Probability density
stairs([xout - dx / 2, xout(end) + dx / 2], [f, f(end)])

xlabel('X'), ylabel('probability density')
% title('Observed maximum of n exponentially distributed variables')

[f, xi] = ksdensity(x);
hold all
plot(xi, f)
hold off

%x = linspace(0,max(xout),500);
%for k = 1:length(n)
%    line(x,gevpdf(x,0,1,log(n(k))),'color','k')
% Matlab <7: use   n(k)*exp(-x).*exp(-n(k)*exp(-x))
%end
