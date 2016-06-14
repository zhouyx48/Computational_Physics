module type_def

type :: vector

real(kind=8) :: x
real(kind=8) :: y

end type vector

end module


program cannon_traj

use type_def

implicit none

! The maximum steps allowed to calculate is 9999.
integer(kind=8),parameter :: max_siz = 9999

! define position vector and velocity vector
type(vector) :: posi(max_siz),velo(max_siz)

! Ohter constsnts (dimension: meter, second)
real(kind=8),parameter :: B_div_m = 4E-05, y_cha = 1E+04, &
& grav_const = 9.8

! time interval delta_t=0.01s
real(kind=8),parameter :: delta_t=0.01

! initial speed v=700m/s
real(kind=8),parameter :: ini_speed = 700., Pi = 3.14159

integer(kind=8) :: aa,bb
real(kind=8) :: time(max_siz)

! initialization, the cannon shell is launched at (0,0), with initial speed v=700m/s, angle=pi/4
time(1) = 0.
posi(1)%x = 0.
posi(1)%y = 0.
velo(1)%x = ini_speed * cos(Pi/4.)
velo(1)%y = ini_speed * sin(Pi/4.)

! Solving ODE
aa = 0
do while (aa < 10000)

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

write(*,*)"Total computation steps:",aa+1

open(10, file="cannon_shell_traj_data_45", form="formatted", &
& access="sequential", recl=493, status="replace")

do bb=1,aa+1

write(10,"(1P5E11.4)")time(bb),posi(bb),velo(bb)

end do

close(10)

end program
