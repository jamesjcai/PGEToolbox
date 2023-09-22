function heateqex2

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


%Solves a sample Dirichlet problem for the heat equation in a rod,

%this time with variable conductivity, 21 mesh points

m = 0; %This simply means geometry is linear.

x = linspace(-5, 5, 21);

t = linspace(0, 4, 81);

sol = pdepe(m, @pdex, @pdexic, @pdexbc, x, t);

% Extract the first solution component as u.

u = sol(:, :, 1);

% A surface plot is often a good way to study a solution.

surf(x, t, u)

title('Numerical solution computed with 21 mesh points in x.')

xlabel('x'), ylabel('t'), zlabel('u')

% A solution profile can also be illuminating.

figure

plot(x, u(end, :))

title('Solution at t = 4')

xlabel('x'), ylabel('u(x,4)')

% --------------------------------------------------------------

    function [c, f, s] = pdex(x, t, u, DuDx)

        c = 1;

        f = (1 + (x / 5).^2) * DuDx; % flux is variable conductivity times u_x

        s = 0;

        % --------------------------------------------------------------

            function u0 = pdexic(x)

                % initial condition at t = 0

                u0 = 20 + 5 * sign(x);

                % --------------------------------------------------------------

                    function [pl, ql, pr, qr] = pdexbc(xl, ul, xr, ur, t)
                        % q's are zero since we have Dirichlet conditions
                        % pl = 0 at the left, pr = 0 at the right endpoint
                        pl = ul - 15;
                        ql = 0;
                        pr = ur - 25;
                        qr = 0;
