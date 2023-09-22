function flignerkilleen(x, g)
% Fligner-Killeen Test of Homogeneity of Variances
%Description
%
%Performs a Fligner-Killeen (median) test of the null that the variances in each of the groups (samples) are the same.
%Usage
%
%fligner.test(x, ...)
%
%## Default S3 method:
%fligner.test(x, g, ...)
%
%## S3 method for class 'formula':
%fligner.test(formula, data, subset, na.action, ...)
%
%Arguments
%x 	a numeric vector of data values, or a list of numeric data vectors.
%g 	a vector or factor object giving the group for the corresponding elements of x. Ignored if x is a list.
%formula 	a formula of the form lhs ~ rhs where lhs gives the data values and rhs the corresponding groups.
%data 	an optional matrix or data frame (or similar: see model.frame) containing the variables in the formula formula. By default the variables are taken from environment(formula).
%subset 	an optional vector specifying a subset of observations to be used.
%na.action 	a function which indicates what should happen when the data contain NAs. Defaults to getOption("na.action").
%... 	further arguments to be passed to or from methods.
%Details
%
%If x is a list, its elements are taken as the samples to be compared for homogeneity of variances, and hence have to be numeric data vectors. In this case, g is ignored, and one can simply use fligner.test(x) to perform the test. If the samples are not yet contained in a list, use fligner.test(list(x, ...)).
%
%Otherwise, x must be a numeric data vector, and g must be a vector or factor object of the same length as x giving the group for the corresponding elements of x.
%
%The Fligner-Killeen (median) test has been determined in a simulation study as one of the many tests for homogeneity of variances which is most robust against departures from normality, see Conover, Johnson & Johnson (1981). It is a k-sample simple linear rank which uses the ranks of the absolute values of the centered samples and weights a(i) = qnorm((1 + i/(n+1))/2). The version implemented here uses median centering in each of the samples (F-K:med X^2 in the reference).
%Value
%
%A list of class "htest" containing the following components:
%statistic 	the Fligner-Killeen:med X^2 test statistic.
%parameter 	the degrees of freedom of the approximate chi-squared distribution of the test statistic.
%p.value 	the p-value of the test.
%method 	the character string "Fligner-Killeen test of homogeneity of variances".
%data.name 	a character string giving the names of the data.
%References
%
%William J. Conover & Mark E. Johnson & Myrle M. Johnson (1981). A comparative study of tests for homogeneity of variances, with applications to the outer continental shelf bidding data. Technometrics 23, 351–361.
%See Also
%
%ansari.test and mood.test for rank-based two-sample test for a difference in scale parameters; var.test and bartlett.test for parametric tests for the homogeneity of variances.
%Examples
%
%require(graphics)
%
%plot(count ~ spray, data = InsectSprays)
%fligner.test(InsectSprays$count, InsectSprays$spray)
%fligner.test(count ~ spray, data = InsectSprays)
%## Compare this to bartlett.test()
%
%
%http://stat.ethz.ch/R-manual/R-patched/library/stats/html/fligner.test.html
%
