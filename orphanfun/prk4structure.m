function [p]=prk4structure(lnp)
%[p]=prk4structure(lnp)
%
%Compute posterior probabilities of K. 
% lnp - ln Pr(X|K), e.g., lnp=[-4356 -3983 -3982 -3983 -4006];
% p - posterior probabilities of K
%
% Assuming a uniform prior on K = {1,...,5}. Then from Bayes' Rule, Pr(K =
% 2) is given by ...

lnp2=lnp-max(lnp);
p=exp(lnp2)./sum(exp(lnp2));


