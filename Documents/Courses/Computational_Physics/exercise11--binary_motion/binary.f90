module type_def

type :: vector

real(kind=8) :: x
real(kind=8) :: y

end type vector

end module

!----------------------------!

program binary

use type_def

implicit none

integer(kind=8) :: aa,cc

! The maximum steps allowed to calculate is 5000.
integer(kind=8),parameter :: max_siz = 5000

! Claim constsnts. Time interval delta_t=0.001 year, star1 have mass M_sun.
real(kind=8),parameter :: Pi = 3.141592654, delta_t = 0.001, m_star1 = 1.0D+00

real(kind=8) :: time(max_siz),r(max_siz),m_star2(4)

! define position vector and velocity vector of two stars
type(vector) :: posi1(max_siz),velo1(max_siz), &
& posi2(max_siz),velo2(max_siz)

m_star2 = (/0.25D+00, 0.5D+00, 0.75D+00, 1.0D+00/)

open(10, file="binary_motion_data", form="formatted", &
& access="sequential", recl=493, status="replace")

do cc = 1,4

write(10,"(1P2E11.4)")m_star1,m_star2(cc)

! initialization.
! star1 is located at (-0.5 AU,0), with velocity (0,-Pi AU/yr).
! star2 is located at (0.5 AU,0), with velocity (0,Pi AU/yr).
posi1(1)%x = -0.5D+00
posi1(1)%y = 0.0D+00
posi2(1)%x = 0.5D+00
posi2(1)%y = 0.0D+00
velo1(1)%x = 0.0D+00
velo1(1)%y = -Pi
velo2(1)%x = 0.0D+00
velo2(1)%y = Pi
time(1) = 0.0D+00
r(1) = sqrt((posi2(1)%y - posi1(1)%y)**2 + &
& (posi2(1)%x - posi1(1)%x)**2)

! Solving ODE
  do aa = 1,max_siz

  velo1(aa+1)%x = velo1(aa)%x - (4. * Pi**2 * m_star2(cc))/(r(aa))**3 * &
& (posi1(aa)%x - posi2(aa)%x) * delta_t
  posi1(aa+1)%x = posi1(aa)%x + velo1(aa+1)%x * delta_t
  velo1(aa+1)%y = velo1(aa)%y - (4. * Pi**2 * m_star2(cc))/(r(aa))**3 * &
& (posi1(aa)%y - posi2(aa)%y) * delta_t
  posi1(aa+1)%y = posi1(aa)%y + velo1(aa+1)%y * delta_t

  velo2(aa+1)%x = velo2(aa)%x - (4. * Pi**2 * m_star1)/(r(aa))**3 * &
& (posi2(aa)%x - posi1(aa)%x) * delta_t
  posi2(aa+1)%x = posi2(aa)%x + velo2(aa+1)%x * delta_t
  velo2(aa+1)%y = velo2(aa)%y - (4. * Pi**2 * m_star1)/(r(aa))**3 * &
& (posi2(aa)%y - posi1(aa)%y) * delta_t
  posi2(aa+1)%y = posi2(aa)%y + velo2(aa+1)%y * delta_t

  r(aa+1) = sqrt((posi2(aa+1)%y - posi1(aa+1)%y)**2 + &
& (posi2(aa+1)%x - posi1(aa+1)%x)**2)

  time(aa+1) = time(aa) + delta_t

  write(10,"(1P6E11.4)")time(aa),r(aa),posi1(aa)%x,posi1(aa)%y, &
& posi2(aa)%x,posi2(aa)%y

  end do

end do

close(10)

end program
