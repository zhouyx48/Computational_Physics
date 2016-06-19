program three_body

use type_def

implicit none

integer(kind=8) :: aa,cc

! The maximum steps allowed to calculate is 5000.
integer(kind=8),parameter :: max_siz = 5000

! Claim constsnts. Time interval delta_t=0.001 year, star1 have mass M_sun.
real(kind=8),parameter :: Pi = 3.141592654, delta_t = 0.001, &
& m_starS = 1.0D+00, m_starE = 3.003D-06

real(kind=8) :: r_EJ(max_siz),r_ES(max_siz),r_JS(max_siz), &
& time(max_siz),m_starJ(4)

! define position vector and velocity vector of two stars
type(vector) :: posiE(max_siz),veloE(max_siz), &
& posiJ(max_siz),veloJ(max_siz),posiS(max_siz)

m_starJ = (/1./1047., 10./1047., 100./1047., 1000./1047./)

open(10, file="three_body_data", form="formatted", &
& access="sequential", recl=493, status="replace")

do cc = 1,4

write(10,"(1P3E11.4)")m_starE,m_starJ(cc),m_starS

! initialization.
! velocity of Earth - (0,2Pi AU/yr); velocity of Jupiter - (0,2.755 AU/yr).
posiE(1)%y = 0.0D+00
posiJ(1)%y = 0.0D+00
veloE(1)%x = 0.0D+00
veloE(1)%y = 2.0D+00 * Pi
veloJ(1)%x = 0.0D+00
veloJ(1)%y = 2.755D+00

! center of mass is always origin point (0,0)
posiS(1)%x = -(m_starE + 5.2D+00 * m_starJ(cc))/ &
& (m_starS + m_starE + m_starJ(cc))
posiS(1)%y = (m_starE * posiE(1)%y + m_starJ(cc) * posiJ(1)%y)/m_starS

! Sun-Earth distance is 1AU while Sun-Jupiter is 5.2AU
posiE(1)%x = posiS(1)%x + 1.0D+00
posiJ(1)%x = posiS(1)%x + 5.2D+00

r_EJ(1) = sqrt((posiJ(1)%y - posiE(1)%y)**2 + &
& (posiJ(1)%x - posiE(1)%x)**2)
r_ES(1) = sqrt((posiS(1)%y - posiE(1)%y)**2 + &
& (posiS(1)%x - posiE(1)%x)**2)
r_JS(1) = sqrt((posiS(1)%y - posiJ(1)%y)**2 + &
& (posiS(1)%x - posiJ(1)%x)**2)

time(1) = 0.0D+00

  ! Solving ODE
  do aa = 1,max_siz

  ! Earth
  veloE(aa+1)%x = veloE(aa)%x - &
& (4. * Pi**2 * m_starJ(cc))/(r_EJ(aa))**3 * &
& (posiE(aa)%x - posiJ(aa)%x) * delta_t - &
& (4. * Pi**2 * m_starS)/(r_ES(aa))**3 * &
& (posiE(aa)%x - posiS(aa)%x) * delta_t

  posiE(aa+1)%x = posiE(aa)%x + veloE(aa+1)%x * delta_t

  veloE(aa+1)%y = veloE(aa)%y - &
& (4. * Pi**2 * m_starJ(cc))/(r_EJ(aa))**3 * &
& (posiE(aa)%y - posiJ(aa)%y) * delta_t - &
& (4. * Pi**2 * m_starS)/(r_ES(aa))**3 * &
& (posiE(aa)%y - posiS(aa)%y) * delta_t

  posiE(aa+1)%y = posiE(aa)%y + veloE(aa+1)%y * delta_t

  ! Jupiter
  veloJ(aa+1)%x = veloJ(aa)%x - &
& (4. * Pi**2 * m_starE)/(r_EJ(aa))**3 * &
& (posiJ(aa)%x - posiE(aa)%x) * delta_t - &
& (4. * Pi**2 * m_starS)/(r_JS(aa))**3 * &
& (posiJ(aa)%x - posiS(aa)%x) * delta_t

  posiJ(aa+1)%x = posiJ(aa)%x + veloJ(aa+1)%x * delta_t

  veloJ(aa+1)%y = veloJ(aa)%y - &
& (4. * Pi**2 * m_starE)/(r_EJ(aa))**3 * &
& (posiJ(aa)%y - posiE(aa)%y) * delta_t - &
& (4. * Pi**2 * m_starS)/(r_JS(aa))**3 * &
& (posiJ(aa)%y - posiS(aa)%y) * delta_t

  posiJ(aa+1)%y = posiJ(aa)%y + veloJ(aa+1)%y * delta_t

  ! coordinate of the Sun
  posiS(aa+1)%x = -(m_starE * posiE(aa+1)%x + &
& m_starJ(cc) * posiJ(aa+1)%x)/m_starS
  posiS(aa+1)%y = -(m_starE * posiE(aa+1)%y + &
& m_starJ(cc) * posiJ(aa+1)%y)/m_starS

  ! distance between the three body
  r_EJ(aa+1) = sqrt((posiJ(aa+1)%y - posiE(aa+1)%y)**2 + &
& (posiJ(aa+1)%x - posiE(aa+1)%x)**2)
  r_ES(aa+1) = sqrt((posiS(aa+1)%y - posiE(aa+1)%y)**2 + &
& (posiS(aa+1)%x - posiE(aa+1)%x)**2)
  r_JS(aa+1) = sqrt((posiS(aa+1)%y - posiJ(aa+1)%y)**2 + &
& (posiS(aa+1)%x - posiJ(aa+1)%x)**2)

  time(aa+1) = time(aa) + delta_t

  write(10,"(1P6E11.4)")time(aa),r_EJ(aa),r_ES(aa),r_JS(aa),posiE(aa)%x,posiE(aa)%y

  end do

end do

close(10)

end program
