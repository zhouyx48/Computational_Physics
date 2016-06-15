module user_supplied_ODE

save

! the value of S_0/m is chosen from textbook p38
! NO SPIN => angular velosity omega = 0rpm
! time interval delta_t = 0.01s
real(kind=8),parameter :: S_div_m = 4.1E-04, grav_const = 9.8, &
& omega = 0., delta_t = 0.01

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

!  Notice that there is no force act on z-direction baseball
!move freely along z-axis. (six ODEs reduce to four.)
real ( kind = 8 ) t
real ( kind = 8 ) y(4)
real ( kind = 8 ) yp(4)

! Following eq. 2.26, coefficient 
! B_2/m = 0.0039 + 0.0058/(1.+exp((v-35.)/5.))
yp(1) = y(2)
yp(2) = - (0.0039 + 0.0058/(1. + &
& exp((sqrt(y(2)**2 + y(4)**2)-35.)/5.))) * &
& sqrt(y(2)**2 + y(4)**2) * y(2)
yp(3) = y(4)
yp(4) = - grav_const + S_div_m * y(2) * omega

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
  integer ( kind = 4 ), parameter :: neqn = 4

  real ( kind = 8 ) abserr
  integer ( kind = 4 ) flag
! integer ( kind = 4 ) i_step
! integer ( kind = 4 ) n_step
! external r8_f1
! real ( kind = 8 ) r8_y1x
  real ( kind = 8 ) relerr
  real ( kind = 8 ) t
  real ( kind = 8 ) t_out
! real ( kind = 8 ) t_start
! real ( kind = 8 ) t_stop
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  abserr = sqrt ( epsilon ( abserr ) )
  relerr = sqrt ( epsilon ( relerr ) )

  flag = 1

! t_start = 0.0D+00
! t_stop = 20.0D+00

! n_step = 5

! initialization. initial position (x(t=0),y(t=0))=(0,1.5m),
! initial speed (v_x(t=0),v_y(t=0)) = (40[m/s],0)
  t_out = 0.0D+00
  t = t_out
  y(1) = 0.0D+00
  y(2) = 40.0D+00
  y(3) = 1.5D+00
  y(4) = 0.0D+00

  call r8_f1 ( t, y, yp )

!  write ( *, '(a)' ) ' '
!  write ( *, '(a)' ) '  FLAG     T             Y            Y''           Y_Exact         Error'
!  write ( *, '(a)' ) ' '
!  write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1), yp(1), r8_y1x ( t ), &
!  y(1) - r8_y1x ( t )

  open(10, file="trajectory_of_nospin_baseball", form="formatted",&
& access="sequential", recl=493, status="replace")

  do

    t = t_out

    t_out = t + delta_t

!   use the rkf45 algorithm to solve our problem
    call r8_rkf45 ( r8_f1, neqn, y, yp, t, t_out, relerr, abserr, flag )

    write (10,'(i4,2x,5g14.6)') flag, t, y(1:4)

      if (y(3)<=0.) then

      exit

      end if

  end do

  close(10)

  return

end program
