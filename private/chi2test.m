function [P,P2,X2,X2C] = chi2test (a,b,c,d);

x=[a,b;c,d];
e=sum(x')'*sum(x)/sum(sum(x));
X2=(x-e).^2./e;
X2=sum(sum(X2));
df=prod(size(x)-[1,1]);
P=1-chi2cdf(X2,df);


% Yate's correction: When there are only two categories  (e.g. male/female) or, more correctly, when there is only one degree of freedom, the c2 test should not, strictly, be used. There have been various attempts to correct this deficiency, but the simplest is to apply Yates correction to our data. To do this, we simply subtract 0.5 from each calculated value of "O-E", ignoring the sign (plus or minus).
X2C=((abs(x-e)-0.5)).^2./e;
X2C=sum(sum(X2C));
P2=1-chi2cdf(X2C,df);


if (nargout<1)
disp(' ')
disp('Table of the Chi-square test results for:')
disp(['a = ' num2str(a) ', b = ' num2str(b) ', c = ' num2str(c) ', d = ' num2str(d) '']);
fprintf('----------------------------------------------\n');
fprintf(['Chi-square value: %10.3f\n'], X2);
fprintf(['      P-value: %10.5f\n'], P);
fprintf(['Chi-square with Yates'' correction: %10.3f\n'], X2C);
fprintf(['      P-value: %10.5f\n'], P2);
fprintf('----------------------------------------------\n');
end