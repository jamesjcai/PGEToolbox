function [y, z, sol] = kimura62(x, gamma, h)
% Fixation probability by solving Kol-backward equation
%
% Kimura 62
% http://www.genetics.org/cgi/reprint/47/6/713
%
% Watterson 1995 - equ (14)
%  Motoo Kimura's use of diffusion theory in pop gen

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin < 3
    h = 0.5;
end


[sol] = ode45(@vdp2, [0, 1], [0; 1], [], gamma, h);
Sxint1 = deval(sol, 1);
Sxint = deval(sol, x);
y = Sxint(1, :) ./ Sxint1(1, 1);
z = Sxint(2, :) ./ Sxint1(2, 1);

%{
[t,y] = ode45(@vdp2,[0 1],[0; 1],[],gamma,h);
plot(t,y(:,1),'-',t,y(:,2),'--')
title('Solution');
xlabel('frequency p');
ylabel('solution u');
legend('u_1','u_2')
vline(x);
hline(Sxint(1,:));
%}

    function dydt = vdp1(t, y, gamma, h)
        dydt = [y(2); -2 * gamma * y(2)];

            function dydt = vdp2(t, y, gamma, h)
                dydt = [y(2); -4 * gamma * (h + (1 - 2 * h) * t) * y(2)];
