function [P] = genedrift(npop,p0,nrepl,ngenr)
%GENEDRIFT - Simulation of genetic drift
% to illustrate the consequences of genetic drift

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if (nargin<1), npop=100; end    % Population size
if (nargin<2), p0=0.5; end      % Initial frequency of A1 allele
if (nargin<3), nrepl=5; end     % Number of replicates
if (nargin<4), ngenr=100; end   % Number of generation

P=zeros(nrepl,ngenr);

for i=1:nrepl
p=p0;
	for j=1:ngenr
		N=2*npop;
		p=sum(rand(1,N)<p)/N;
		P(i,j)=p;
	end
end

if (nargout<1),
    subplot(2,2,1)
hold on
      for k=1:nrepl, plot(1:ngenr,P(k,:)); end
      xlabel('Generations')
      ylabel('Allele Frequency p')
      title(sprintf('Change of allele frequency for population (diploid) of size N = %d', npop))
      ylim([0 1])
hold off
subplot(2,2,3)
H=2*P.*(1-P);
hold on
for k=1:nrepl, plot(1:ngenr,H(k,:),'r'); end
      xlabel('Generations')
      ylabel('Het')
      title(sprintf('Change of heterozygosity of population'))
      ylim([0 0.52])
hold off

end


 
 

%drift_markov(p0,npop,ngenr);
% drift_markov(p_start,N,generations)
%
% This function simulates genetic drift as a Markov process.  This is to say, a transition matrix is
% created where the Pr(i,j) defines the probability of change in a generation from the state i
% to the state j.  This transition matrix is created and the the population shifts between states
% according to their probability defined by the matrix.  This is not remarkable for producing a
% different outcome than genetic drift by sampling, but for producing the same outcome.
% Input initial p (p_start), population size (N), and number of generations (generations).
% Note that the population is a population of N haploids.
%
% Function written by Liam Revell 2005 for Matlab R12.


function drift_markov(p_start,N,generations) %#ok<*DEFNU>

% create the transition matrix
for i=0:N
    for j=0:N
        Pr(i+1,j+1)=(1/(factorial(N-i)*factorial(i)))*(j/N)^i*(1-j/N)^(N-i)*(factorial(N));
    end
end

time(1)=0;
Ni=double(int32(p_start*N));
for i=2:generations
    test=0;
    while test==0
        % pick a new random frequency between 0 and N
        newNi=int32(rand*(N+1));
        % change to that state with the probability defined by Pr[N(i),N(i-1)]
        if Pr(double(newNi)+1,double(Ni(i-1))+1)>rand
            Ni(i)=newNi;
            test=1;
        end
    end
    time(i)=i-1;
end
p=Ni/N;

% plot p over time
plot(time,p,'k');
hold on;
axis([0 generations 0 1]);
xlabel('Time');
ylabel('p_A');
title('Drift as a Markov process');





% drift_variance(p_start,N,generations)
%
% This function plots the expected change in variance in p accross populations over time.
% It also creates an animated histogram of allele number by probability.  This produces a neat
% animation for small haploid N (say N=10-40, generations 20-60).
% Input initial p (p_start), population size (N), and number of generations (generations).
% Note that the population is a population of N haploids.
% Compare histograms to Fig 3.1 in Rice 2004.
%
% Function written by Liam Revell 2005 for Matlab R12.

function drift_variance(p_start,N,generations)

% Calculate transition matrix
for i=0:N
    for j=0:N
        Pr(i+1,j+1)=(1/(factorial(N-i)*factorial(i)))*(j/N)^i*(1-j/N)^(N-i)*(factorial(N));
    end
    N_A(i+1)=0;
    if(p_start==(i/N))
        N_A(i+1)=1;
    end
    var_vect(i+1)=(i/N-p_start)^2;
end
var(1)=0;
time(1)=0;

% Calculate variance over time and expected allele frequencies
for i=1:generations
    N_At(i,:)=N_A;
    time(i+1)=i;
    newN=N_A*Pr';
    var(i+1)=N_A*Pr'*var_vect';
    N_A=newN;
end


plot(time,var);
hold on;
axis([0 generations, 0 0.3]);
xlabel('Generations');
ylabel('var(p_A)');
hold off;

figure;
x=0:1:N;
for i=1:generations
    h=bar(x,N_At(i,:));
    hold on;
    axis([0 N,0 1]);
    xlabel('Number of A alleles');
    ylabel('Frequency');
    hold off;
    M(i)=getframe;
end
movie(M);
