subroutine p11_f ( m, n, x, f )
!
!*******************************************************************************
!
!! P11_F evaluates the M nonlinear functions for problem 11.
!
!
!  Modified:
!
!    28 October 2000
!
!  Parameters:
!
!    Input, integer M, the number of functions.
!
!    Input, integer N, the number of unknowns.
!
!    Input, real X(N), the point at which F is to be evaluated.
!
!    Output, real F(M), the value of the functions evaluated at X.
!
  implicit none
!
  integer m
  integer n
!
  real div
  real dx
  real f(m)
  integer i
  integer j
  real s1
  real s2
  real x(n)
!
  do i = 1, 29

    div = real ( i ) / 29.0E+00
    s1 = 0.0E+00
    dx = 1.0E+00
    do j = 2, n
      s1 = s1 + real ( j - 1 ) * dx * x(j)
      dx = div * dx
    end do

    s2 = 0.0E+00
    dx = 1.0E+00
    do j = 1, n
      s2 = s2 + dx * x(j)
      dx = div * dx
    end do
    f(i) = s1 - s2**2 - 1.0E+00
  end do

  f(30) = x(1)
  f(31) = x(2) - x(1)**2 - 1.0E+00

  return
end
subroutine p11_j ( ldfjac, m, n, x, fjac )
!
!*******************************************************************************
!
!! P11_J evaluates the jacobian for problem 11.
!
!
!  Modified:
!
!    28 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDFJAC, the leading dimension of FJAC.
!
!    Input, integer M, the number of equations.
!
!    Input, integer N, the number of variables.
!
!    Input, real X(N), the point at which the jacobian is to be evaluated.
!
!    Output, real FJAC(LDFJAC,N), the M by N jacobian matrix.
!
  implicit none
!
  integer ldfjac
  integer n
!
  real div
  real dx
  real fjac(ldfjac,n)
  integer i
  integer j
  integer m
  real s2
  real temp
  real v(11)
  real x(n)
!
  do i = 1, 29
    div = real ( i ) / 29.0E+00
    s2 = 0.0E+00
    dx = 1.0E+00
    do j = 1, n
      s2 = s2 + dx * x(j)
      dx = div * dx
    end do
    temp = 2.0E+00 * div * s2
    dx = 1.0E+00 / div
    do j = 1, n
      fjac(i,j) = dx * ( real ( j - 1 ) - temp )
      dx = div * dx
    end do
  end do

  fjac(30:31,1:n) = 0.0E+00
  fjac(30,1) = 1.0E+00
  fjac(31,1) = -2.0E+00 * x(1)
  fjac(31,2) = 1.0E+00

  return
end
subroutine p11_sol ( m, n, known, x )
!
!*******************************************************************************
!
!! P11_SOL returns the solution of problem 11.
!
!
!  Modified:
!
!    28 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of equations.
!
!    Input, integer N, the number of variables.
!
!    Output, integer KNOWN, 1 or 0, if the solution is known or not.
!
!    Output, real X(N), the solution, if known.
!
  implicit none
!
  integer n
!
  integer known
  integer m
  real x(n)
!
  known = 0
  x(1:n) = 0.0E+00

  return
end
subroutine p11_start ( n, x )
!
!*******************************************************************************
!
!! P11_START sets a starting point for problem 11.
!
!
!  Modified:
!
!    28 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of variables.
!
!    Output, real X(N), a starting point for the problem.
!
  implicit none
!
  integer n
!
  real x(n)
!
  x(1:n) = 0.0E+00

  return
end
subroutine p11_title ( title )
!
!*******************************************************************************
!
!! P11_TITLE specifies the title for problem 11.
!
!
!  Modified:
!
!    28 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the problem title.
!
  character ( len = * ) title
!
  title = '11: Watson function.'

  return
end
