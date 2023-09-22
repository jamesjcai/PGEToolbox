function brunnermunzel(x, y, alternative, alpha)
%brunnermunzel - Brunner-Munzel Test for Stochastic Equality
%
%x 	the numeric vector of data values from the sample 1.
%y 	the numeric vector of data values from the sample 2.
%alpha 	confidence level, default is 0.05 for 95 interval.
%alternative 	a character string specifying the alternative hypothesis, must be one of 'two.sided' (default), 'greater' or 'less'. User can specify just the initial letter.
%
%statistic 	the Brunner-Munzel test statistic.
%parameter 	the degrees of freedom.
%conf.int 	the confidence interval.
%p.value 	the p-value of the test.
%data.name 	a character string giving the name of the data.
%
%References:
%Brunner, E. and Munzel, U. (2000) The Nonparametric Behrens-Fisher Problem: Asymptotic Theory and a Small-Sample Approximation, Biometrical Journal 42, 17-25.
%Reiczigel, J., Zakarias, I. and Rozsa, L. (2005) A Bootstrap Test of Stochastic Equality of Two Populations, The American Statistician 59, 1-6.
%
%See Also: RANKSUM, KSTEST2
%
%Examples
%
%## Pain score on the third day after surgery for 14 patients under
%## the treatment \emph{Y} and 11 patients under the treatment \emph{N}
%## (see Brunner and Munzel (2000))
%
%Y<-c(1,2,1,1,1,1,1,1,1,1,2,4,1,1)
%N<-c(3,3,4,3,1,2,3,1,1,5,4)
%
%brunner.munzel.test(Y, N)
%
%##       Brunner-Munzel Test
%## data: Y and N
%## Brunner-Munzel Test Statistic = 3.1375,  df = 17.683, p-value = 0.005786
%## 95 percent confidence interval:
%##  0.5952169 0.9827052
%## sample estimates:
%
% http://rss.acs.unt.edu/Rdoc/library/lawstat/html/brunner.munzel.test.html
% http://www3.interscience.wiley.com/journal/122219068/abstract
% DOI: 10.1002/sim.3561


## P(X < Y) + .5 * P(X = Y)
## 0.788961
