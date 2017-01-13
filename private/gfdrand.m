function [y]=gfdrand(x,t,x0)
% gfdrand - gene frequency distribution under random genetic drift

r=1-2*x0;
z=1-2*x;

j=1; y=0;
while (j<100)
    A=((2*j+1).*(1-r.^2))./(j.*(j+1));
    B=gegenbauerc(j,1.5,r);
    C=gegenbauerc(j,1.5,z);
    D=exp(-0.5*j.*(j+1).*t);
    incr = A.*B.*C.*D;
    j=j+1;
    y=y+incr;
end

y=4*x0.*(1-x0).*y;




function [c]=gegenbauerc(n,m,x)
%GegenbauerC[n, m, x] gives the n-th Gegenbauer polynomial in x for parameter m.

cx=gegenbauer_poly(n,m,x);
c=cx(end);


function cx = gegenbauer_poly(n,alpha,x)
% GEGENBAUER_POLY computes the Gegenbauer polynomials C(I,ALPHA)(X).
%
%  Discussion:
%
%    The Gegenbauer polynomial can be evaluated in Mathematica with
%    the command
%
%      GegenbauerC[n,m,x]
%
%  Differential equation:
%
%    (1-X*X) Y'' - (2 ALPHA + 1) X Y' + N (N + 2 ALPHA) Y = 0
%
%  Recursion:
%
%    C(0,ALPHA,X) = 1,
%    C(1,ALPHA,X) = 2*ALPHA*X
%    C(N,ALPHA,X) = (  ( 2*N-2+2*ALPHA) * X * C(N-1,ALPHA,X) 
%                    + (  -N+2-2*ALPHA)   *   C(N-2,ALPHA,X) ) / N
%
%  Restrictions:
%
%    ALPHA must be greater than -0.5.
%
%  Special values:
%
%    If ALPHA = 1, the Gegenbauer polynomials reduce to the Chebyshev
%    polynomials of the second kind.
%
%  Norm:
%
%    Integral ( -1 <= X <= 1 ) 
%      ( 1 - X**2 )**( ALPHA - 0.5 ) * C(N,ALPHA,X)**2 dX
%
%    = PI * 2**( 1 - 2 * ALPHA ) * Gamma ( N + 2 * ALPHA ) 
%      / ( N! * ( N + ALPHA ) * ( Gamma ( ALPHA ) )**2 )
%
%  Modified:
%
%    06 August 2004
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Stephen Wolfram,
%    The Mathematica Book,
%    Fourth Edition,
%    Wolfram Media / Cambridge University Press, 1999.
%
%  Parameters:
%
%    Input, integer N, the highest order polynomial to compute.
%    Note that polynomials 0 through N will be computed.
%
%    Input, real ALPHA, a parameter which is part of the definition of
%    the Gegenbauer polynomials.  It must be greater than -0.5.
%
%    Input, real X, the point at which the polynomials are to be evaluated.
%
%    Output, real CX(1:N+1), the values of the first N+1 Gegenbauer
%    polynomials at the point X.  
%
  if ( alpha <= -0.5 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'GEGENBAUER_POLY - Fatal error!\n' );
    fprintf ( 1, '  Illegal value of ALPHA = %f\n', alpha );
    fprintf ( 1, '  but ALPHA must be greater than -0.5.\n' );
    error ( 'GEGENBAUER_POLY - Fatal error!' );
  end

  if ( n < 0 )
    cx = [];
    return
  end

  cx(1) = 1.0;

  if ( n == 0 )
    return
  end

  cx(2) = 2.0E+00 * alpha * x;

  for i = 2 : n
    cx(i+1) = (  (     2 * i - 2  + 2.0 * alpha ) * x * cx(i)     ...
             +   (       - i + 2  - 2.0 * alpha ) *     cx(i-1) ) ...
             /             i ;
  end
  