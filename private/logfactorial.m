function p = logfactorial(n)
%LOGFACTORIAL Log Factorial function.
%   FACTORIAL(N) is the sum of all the integers from 1 to N,
%   i.e. prod(1:N). Since double precision numbers only have about
%   15 digits, the answer is only accurate for N <= 21. For larger N,
%   the answer will have the right magnitude, and is accurate for
%   the first 15 digits.
%
%   See also FACTORIAL.

if (length(n) ~= 1) | (fix(n) ~= n) | (n < 0)
    error('MATLAB:logfactorial:NNotPositiveInteger', ...
        'N must be a positive integer.');
end

p = sum(log(1:n));