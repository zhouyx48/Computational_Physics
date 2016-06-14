module type_def

type :: vector

real(kind=8) :: x
real(kind=8) :: y

end type vector

end module

!----------------------------!

module get_data

use type_def

save

integer(kind=8) :: aa,bb,cc

! The maximum steps allowed to calculate is 49999.
integer(kind=8),parameter :: max_siz = 49999

contains

function dat(theta)

implicit none

real(kind=8) :: time(max_siz),dat(max_siz,5),theta

! define position vector and velocity vector
type(vector) :: posi(max_siz),velo(max_siz)

! Ohter constsnts (dimension: meter, second)
real(kind=8),parameter :: B_div_m = 4E-05, y_cha = 1E+04, &
& grav_const = 9.8

! time interval delta_t=0.01s
real(kind=8),parameter :: delta_t=0.01

! initial speed v=700m/s
real(kind=8),parameter :: ini_speed = 700., Pi = 3.14159

! initialization, the cannon shell is launched at (0,0), with initial speed v=700m/s, angle=pi/4
time(1) = 0.
posi(1)%x = 0.
posi(1)%y = 0.
velo(1)%x = ini_speed * cos((Pi*theta)/180.)
velo(1)%y = ini_speed * sin((Pi*theta)/180.)

! Solving ODE
aa = 0
do while (aa < max_siz+1)

aa = aa+1

posi(aa+1)%x = posi(aa)%x + velo(aa)%x * delta_t
velo(aa+1)%x = velo(aa)%x - B_div_m * exp(-posi(aa)%y/y_cha) *&
& sqrt(velo(aa)%x**2 + velo(aa)%y**2) * velo(aa)%x * delta_t

posi(aa+1)%y = posi(aa)%y + velo(aa)%y * delta_t
velo(aa+1)%y = velo(aa)%y - grav_const * delta_t -&
& B_div_m * exp(-posi(aa)%y/y_cha) * &
& sqrt(velo(aa)%x**2 + velo(aa)%y**2) * velo(aa)%y * delta_t

time(aa+1) = time(aa)+delta_t

  if (posi(aa+1)%y <= 0.) then

  exit

  end if

end do

dat(:,1) = time(:)
dat(:,2) = posi(:)%x
dat(:,3) = posi(:)%y
dat(:,4) = velo(:)%x
dat(:,5) = velo(:)%y

end function

end module get_data

!--------------------!

program find_max_range

use type_def
use get_data

implicit none

! define position vector and velocity vector
type(vector) :: posi,velo

! Ohter constsnts (dimension: meter, second)
real(kind=8),parameter :: B_div_m = 4E-05, y_cha = 1E+04, &
& grav_const = 9.8

! time interval delta_t=0.01s
real(kind=8),parameter :: delta_t=0.01

! initial speed v=700m/s
real(kind=8),parameter :: ini_speed = 700., Pi = 3.14159

real(kind=8) :: max_angle,max_range,ang
real(kind=8) :: result(max_siz,5)

max_range = 0.
max_angle = 0.

! scan to find the launching angle corresponds to maximum range.
do cc=1,899

ang = real(cc)/10.

! initialization, the cannon shell is launched at (0,0), 
! with initial speed v=700m/s, angle=pi/4
posi%x = 0.
posi%y = 0.
velo%x = ini_speed * cos((Pi/180.)*ang)
velo%y = ini_speed * sin((Pi/180.)*ang)

  ! Solving ODE
  do

  posi%x = posi%x + velo%x * delta_t
  velo%x = velo%x - B_div_m * exp(-posi%y/y_cha) *&
& sqrt(velo%x**2 + velo%y**2) * velo%x * delta_t

  posi%y = posi%y + velo%y * delta_t
  velo%y = velo%y - grav_const * delta_t -&
& B_div_m * exp(-posi%y/y_cha) * &
& sqrt(velo%x**2 + velo%y**2) * velo%y * delta_t

    if (posi%y <= 0.) then

      if (posi%x >= max_range) then

      max_range = posi%x
      max_angle = ang

      end if

    exit

    end if

  end do

end do

write(*,*)"Maximum range of the cannon is",max_range, &
& "coreresponds to theta=",max_angle

result = dat(max_angle)

open(10, file="cannon_shell_traj_data_max_range", form="formatted", &
& access="sequential", recl=493, status="replace")

do bb=1,aa+1

write(10,"(1P5E11.4)")result(bb,1:5)

end do

close(10)

end program
