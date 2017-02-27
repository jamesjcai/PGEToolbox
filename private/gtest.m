function [P,G] = gtest (a,b,c,d);

  x=[a,b;c,d];
  [r,c] = size(x);
  if (min([r,c])<2)
    error('Error: two-dimensional table needed for test');
  end;

  R = sum(x);
  C = sum(x')';
  N = sum(R);

  f = x(:);
  G = 2*(sum(f.*log(f)) - sum(R.*log(R)) - sum(C.*log(C)) + N*log(N));
  df = (r-1)*(c-1);

  P = 1-chi2cdf(G,df);


if (nargout<1)
disp(' ')
disp('Table of the G test results for:')
disp(['a = ' num2str(a) ', b = ' num2str(b) ', c = ' num2str(c) ', d = ' num2str(d) '']);
fprintf('----------------------------------------------\n');
fprintf(['G value: %10.3f\n'], G);
fprintf(['   P-value: %10.5f\n'], P);
fprintf('----------------------------------------------\n');
end