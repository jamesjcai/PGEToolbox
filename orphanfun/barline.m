function barline(p, pos)
%barline(p,pos)
%

if nargin < 2
    pos = 1:length(p);
end


%evenpos=min(pos)+[0:length(pos)-1].*((max(pos)-min(pos))/(length(pos)-1));

hold on
nb = length(p);
for ii = 1:nb
    xx = [pos(ii), pos(ii)];
    %xx = [evenpos(ii) pos(ii)];
    yy = [0, p(ii)];
    plot(xx, yy, '-')
end

%ylabel('DAF')
%xlabel('Position on Chromosome')
%plot([min(min(pos)) max(max(pos))],[1 1],'-r')
