function [u]=fixprob(p,s,Ne)
%FIXPROB - fixation probability of a mutation
%
%Usage: [u]=fixprob(p,s,N)
% http://www.genetics.org/cgi/reprint/47/6/713
%
% Kimura 1962 - equ (8)
% http://www.genetics.org/cgi/reprint/47/6/713

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

%{
figure;
N=10000;
[s,p] = meshgrid(0:0.05:1);
s=5*s./N;
f=fixprob(p,s,N);
surf(s,p,f)
%colorbar
xlabel('N_e=10000')
ylabel('p')
zlabel('Probability of fixation')
shading interp


N=10000;
[s,p] = meshgrid(0:0.001:1);
s=5*s./N-2.5./N;
%f=fixprob(p,s,N);
f=log(fixprob(p,s,N));
figure;
surf(s,p,f)
shading interp


p=0.2
s=linspace(-2,2);
u=fixprob(10000,s./10000,p);
plot(s,u)
p=0.01;
hold on
u=fixprob(10000,s./10000,p);
plot(s,u,'r')
p=0.1;
u=fixprob(10000,s./10000,p);
plot(s,u,'k')


figure; hold all
s=linspace(-0.005,0.005,100);
N=200; plot(s,fixprob(1/(2*N),s,N))
N=600; plot(s,fixprob(1/(2*N),s,N))
N=1000; plot(s,fixprob(1/(2*N),s,N))
ylim([0, 0.01]); xlabel('s')
ylabel('P(fixation of new mutation)')
legend({'N=200','N=600','N=1000'},2)


figure; hold all
s=linspace(-0.005,0.005,50);
N=300; plot(s,2*N*fixprob(1/(2*N),s,N),'-k')
N=200; plot(s,2*N*fixprob(1/(2*N),s,N),'--k')
N=100; plot(s,2*N*fixprob(1/(2*N),s,N),':k')
%ylim([0, 4]); 
xlabel('s')
ylabel('P(fixation of new mutation)*2N')
legend({'N=300','N=200','N=100'},2)
vline(0); hline(1)
%}

if nargin==3, a=2*Ne*s; end
if nargin==2, a=s; end

%if nargin<1, N=10000; end
%if nargin<2, s=0; end
%if nargin<3, p=1./N; end

%p=(1-exp(-f))./(1-exp(-2.*s.*N.*f));

if a==0, u=p; return; end
Q1=1-exp(-2*a.*p);
Q2=1-exp(-2*a);
u=Q1./Q2;

