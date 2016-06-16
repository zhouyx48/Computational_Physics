module user_supplied_ODE

save

! length of the the pendulum is set len = 1m, 
! same as fig. 3.5 in txtbook.
! The amplitude of driving force F_D is also defined as unit.
! time interval delta_t = 0.01s
real(kind=8),parameter :: grav_const = 9.8, len = 1., &
& F_D = 1., delta_t = 0.01

real(kind=8) :: frac_const, Omega_D

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
yp(2) = - (grav_const/len) * y(1) - frac_const * y(2) + &
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

  open(10, file="oscillation_data", form="formatted",&
& access="sequential", recl=493, status="replace")

  do aa = 1,3

  ! determining parameter q
  if (aa==1) then

  frac_const = 0.1D+00

  else if (aa==2) then

  frac_const = 2.0D+00

  else

  frac_const = 10.0D+00

  end if

    ! determind Omega_D
    do bb = 1,3

      if (bb==1) then

      Omega_D = 1.0D+00

      else if (bb==2) then

      Omega_D = 3.130495168

      else

      Omega_D = 8.0D+00

      end if

    write (10,'(2g14.6)')frac_const, Omega_D

    ! initial value theta_0 = 0.2rad is the same as p56 fig. 3.5(a).
    t_out = 0.0D+00
    t = t_out
    y(1) = 0.2D+00
    y(2) = -0.2D+00

    call r8_f1 ( t, y, yp )

      ! 2000 steps.
      do step = 1,2000

      t = t_out

      t_out = t + delta_t

      !   use the rkf45 algorithm to solve our problem
      call r8_rkf45 ( r8_f1, neqn, y, yp, t, t_out, relerr, abserr, flag )

      write (10,'(2i4,2x,3g14.6)')step, flag, t, y(1:neqn)

      end do

    end do

  end do

  close(10)

  return

end program
