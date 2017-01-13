function kimura55
%In the diffusion framework the Kolmogorov
%forward equation (also known as the Fokker–Planck
%equation) describes the evolution ofthe allele frequency
%distribution over time, ... 
%
%This PDE is difficult to solve analytically and to our
%knowledge no analytical solution has been determined
%for the general case of selection, mutation, and drift. In
%this case numerical methods have to be used to study
%the dynamics ofthe allele frequency distribution.
%
%The
%solution at equilibrium is easier to obtain and was given
%by Kimura (1955) confirming the result originally found
%by Wright (1937) using a different approach.
%
%kolm_backward
%
%Kimura, M. 1955. Stochastic processes and distribution of gene frequencies under natural selection, Cold Spring Harbor Symp.Quant.Biol. 20, 33–53.
%
%Wright, S. 1937. The distribution of gene frequencies in populations,
%Genetics 23, 307–320.

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

m = 0;
x = linspace(0,1,20);
t = linspace(1,10,20);


sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);
% Extract the first solution component as u.  This is not necessary
% for a single equation, but makes a point about the form of the output.
u = sol(:,:,1);

% A surface plot is often a good way to study a solution.
figure;
surf(x,t,u); 
title('Numerical solution computed with 20 mesh points.');
xlabel('Frequency x');
ylabel('Time t');


nt = length(t);
I = zeros(1,nt);
J = zeros(1,nt);
iok = 2:nt;
for j = iok
  % At time t(j), compute Du/Dx at x = 0.
  [I(j),J(j)] = pdeval(m,x,u(j,:),x(end-3));  
end
% I(t) = (I_p*d/K)*Du(0,t)/Dx
%I = (I_p*d/K)*I;

%figure;
%plot(t(iok),I(iok));
%plot(x(2:end),J(iok));


figure;
%plot(t(iok),I(iok));
plot(x(2:end),I(iok));
legend('From PDEPE');
title('I(t)');
xlabel('Time t');
hold on
%figure;
plot(x,u(:,end-3),'o');
title('Solutions at f = 0.85.');
legend('Numerical, 20 mesh points',0);
xlabel('time t');
ylabel('u(0.85,t)');



% A solution profile can also be illuminating.
figure;
plot(x,u(end,:),'-or');
hold all
plot(x,u(end-1,:),'-o');
plot(x,u(end-2,:),'-o');
title('Solutions at t = 2.');
legend('Numerical, 20 mesh points',0);
xlabel('frequency x');
ylabel('u(x,2)');


function [c,f,s] = pdex1pde(x,t,u,DuDx)
A=0.5.*x.*(1-x);
B=1-2.*x;
c = 1;
f = A.*DuDx;
s = B.*DuDx-u;

%function [c,f,s] = pdex1pde(x,t,u,DuDx)
%Sx=2;
%V=x.*(1-x);
%M=Sx.*x.*(1-x);
%c = 1;
%f = 0.5*V.*DuDx;
%s = M.*DuDx;


function u0 = pdex1ic(x)
%N=100000;
u0 = x;  % (1/2*N);
    
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = ul-0;
ql = 0;
pr = ur-1;
qr = 0;

