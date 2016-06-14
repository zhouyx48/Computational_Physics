program two_nuclei_decay

implicit none

! decay time constant is set as 1s, time interval
! delta_t=0.01s
real(kind=8),parameter :: time_const_tau=1.,delta_t=0.01
integer(kind=8),parameter :: step=1000
integer(kind=8) :: aa
real(kind=8) :: nuclei_num_A(step),nuclei_num_B(step),time(step)

! initialization, following textbook p17 1.5
nuclei_num_A(1) = 100.
nuclei_num_B(1) = 0.
time(1) = 0.

! Solving ODE
do aa=1,step-1

nuclei_num_A(aa+1) = nuclei_num_A(aa) + &
& (nuclei_num_B(aa)/time_const_tau - &
& nuclei_num_A(aa)/time_const_tau) * delta_t

nuclei_num_B(aa+1) = nuclei_num_B(aa) + &
& (nuclei_num_A(aa)/time_const_tau - &
& nuclei_num_B(aa)/time_const_tau) * delta_t

time(aa+1) = time(aa)+delta_t

end do

open(10, file="two_nuclei_decay_data", form="formatted", &
& access="sequential", recl=493, status="replace")

do aa=1,step-1

write(10,"(1P3E11.4)")time(aa),nuclei_num_A(aa),nuclei_num_B(aa)

end do

close(10)

end program two_nuclei_decay
