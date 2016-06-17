module user_supplied_ODE

save

! All the same as sec. 3.6, sigma = 10 and b = 8/3
! time interval delta_t = 0.0001s
real(kind=8),parameter :: delta_t = 0.0001, sigma = 10.0D+00, b = 8./3.

real(kind=8),parameter :: grav_const = 9.8D+00, Pi = 3.141592654

! ATTENTION: r should be convert to real type when solving ODE.
integer(kind=4) :: r

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
real ( kind = 8 ) y(3)
real ( kind = 8 ) yp(3)

! EOM for Lorenz model
yp(1) = sigma * (y(2) - y(1))
yp(2) = - y(1) * y(3) + real(r, kind=8) * y(1) - y(2)
yp(3) = y(1) * y(2) - b * y(3)

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
  integer ( kind = 4 ), parameter :: neqn = 3

  real ( kind = 8 ) abserr
  integer ( kind = 4 ) flag, step, aa, bb
  real ( kind = 8 ) relerr
  real ( kind = 8 ) t
  real ( kind = 8 ) t_out
  real ( kind = 8 ) y(neqn), pre_y(neqn)
  real ( kind = 8 ) yp(neqn)
  real ( kind = 8 ) y_val_x0, z_val_x0

  abserr = sqrt ( epsilon ( abserr ) )
  relerr = sqrt ( epsilon ( relerr ) )

  flag = 1

  open(10, file="Lorenz_model_data", form="formatted",&
& access="sequential", recl=493, status="replace")
  open(11, file="y-z_phase_diagram(x=0)", form="formatted",&
& access="sequential", recl=493, status="replace")

  do r = 0,160,5

  write(10,'(i3)')r

  ! initial conditions x_0 = 1 and y = z = 0 are the same as p79, fig. 3.16.
  t_out = 0.0D+00
  t = t_out
  y(1) = 1.0D+00
  y(2) = 0.0D+00
  y(3) = 0.0D+00
  pre_y(1) = 1.0D+00
  pre_y(2) = 0.0D+00
  pre_y(3) = 0.0D+00

  call r8_f1 ( t, y, yp )

    ! 300,000 steps (t from 0 by 0.0001 to 30s).
    do step = 1,300000

    t = t_out

    t_out = t + delta_t

    !   use the rkf45 algorithm to solve our problem
    call r8_rkf45 ( r8_f1, neqn, y, yp, t, t_out, relerr, abserr, flag )

      ! collect data for y-z phase diagram (x=0)
      if (pre_y(1) * y(1) <= 0.) then

      y_val_x0 = y(2) - y(1) * (y(2) - pre_y(2))/(y(1) - pre_y(1))
      z_val_x0 = y(3) - y(1) * (y(3) - pre_y(3))/(y(1) - pre_y(1))

      write(11,'(i3,i8,1P3E11.4)')r, step, t, y_val_x0, z_val_x0

      end if

    pre_y = y

      ! there are two many data, we only output a small fraction of them.
      if (mod(step,100) == 0) then

      write (10,'(1i8,1P4E11.4)')step, t, y(1:neqn)

      end if

    end do

  end do

  close(10)
  close(11)

  return

end program
