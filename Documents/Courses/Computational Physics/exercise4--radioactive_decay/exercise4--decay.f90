program decay

implicit none

! decay time constant is set as 1s, time interval
! delta_t=0.01s
real(kind=8),parameter :: time_const_tau=1.,delta_t=0.01
integer(kind=8),parameter :: step=1000
integer(kind=8) :: aa
real(kind=8) :: nuclei_num(step),time(step)

! initialization
nuclei_num(1) = 2000.
time(1) = 0.

! Solving ODE
do aa=1,step-1

nuclei_num(aa+1) = nuclei_num(aa)-&
& (nuclei_num(aa)/time_const_tau)*delta_t

time(aa+1) = time(aa)+delta_t

end do

open(10, file="decay_data", form="formatted", &
& access="sequential", recl=493, status="replace")

do aa=1,step-1

write(10,"(1P2E11.4)")nuclei_num(aa),time(aa)

end do

close(10)

end program decay
