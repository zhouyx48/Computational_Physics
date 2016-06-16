module user_supplied_ODE

save

! All the same as fig. 3.6, 
! length of the the pendulum is set len = 9.8m, 
! friction constant q = 0.5, forcing period Omega_D = 2./3.
! The amplitude of driving force F_D = 1.2, as depiected fig. 3.9.
! time interval delta_t = 0.04s
real(kind=8),parameter :: grav_const = 9.8D+00, len = 9.8D+00, &
& F_D = 1.2D+00, Omega_D = 2./3., fric_const = 0.5D+00, &
& delta_t = 0.04, Pi = 3.141592654

contains

subroutine r8_f1 ( t, y, yp )

!  R8_F1 evaluates the derivative for the ODE.
!
!  The original version was written on 26 March 2004,
!by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the value of the independent variable.
!
!    Input, real ( kind = 8 ) Y(NEQN), the value of the dependent variable.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative
!    dY(1:NEQN)/dT.
!
implicit none

real ( kind = 8 ) t
real ( kind = 8 ) y(2)
real ( kind = 8 ) yp(2)

! EOM for oscillating pendulum
yp(1) = y(2)
yp(2) = - (grav_const/len) * sin(y(1)) - fric_const * y(2) + &
& F_D * sin(Omega_D * t)

return

end subroutine r8_f1

end module


!--------------------------!


program rkf45_solve_ODE
!*****************************************************************************80
!
!  rkf45_solve_ODE solves a scalar ODE system using 'rkf45.f90'.
!
!  The original version was written on 17 June 2006,
!by John Burkardt.

  use user_supplied_ODE

  implicit none

! how many ODEs?
  integer ( kind = 4 ), parameter :: neqn = 2

  real ( kind = 8 ) abserr
  integer ( kind = 4 ) flag, step, aa, bb
  real ( kind = 8 ) relerr
  real ( kind = 8 ) t
  real ( kind = 8 ) t_out
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  abserr = sqrt ( epsilon ( abserr ) )
  relerr = sqrt ( epsilon ( relerr ) )

  flag = 1

  open(10, file="chaos_data", form="formatted",&
& access="sequential", recl=493, status="replace")

  ! initial value theta_0 = 0.2rad and omega_0 = 0 is the same as p63.
  t_out = 0.0D+00
  t = t_out
  y(1) = 0.2D+00
  y(2) = 0.0D+00

  call r8_f1 ( t, y, yp )

    ! 99999 steps.
    do step = 1,99999

    t = t_out

    t_out = t + delta_t

    !   use the rkf45 algorithm to solve our problem
    call r8_rkf45 ( r8_f1, neqn, y, yp, t, t_out, relerr, abserr, flag )

      ! if theta exceed Pi/-Pi, then subtract/add 2Pi
      if (y(1) > Pi) then

      y(1) = y(1) - 2. * Pi

      else if (y(1) < -Pi) then

      y(1) = y(1) + 2. * Pi

      end if

    write (10,'(2i5,2x,3g14.6)')step, flag, t, y(1:neqn)

    end do

  close(10)

  return

end program
