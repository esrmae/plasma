Program cfd_incomp
!
Use shared_data
Use setup 
Use flow_control
Use flow_field
Use boundary_conditions
Use cal_coefficient
Use solvers
Use poisson_solver
Use update_field
Use implicit_solver

Implicit None
Integer :: i,j,k
Double Precision :: start_time, end_time, tot_time
!
Call initia
If(ibcjet == 1) Call initial_control
If (ifdt == 1) Call cal_dt

time0 = time

start_time = MPI_Wtime()
!
Do ii=1,nstep
  If (iconbc == 0 .And. imass == 1) Call compute_dpdx_mean
  If (ibcjet == 1) Call compute_flow_control
!
!   Call compute_boundary
!
!	Implicit solve 
!
    Call implicit_solve
!
!	Poisson equation solution
!
    Call compute_poisson
!
  If (((MPI_Wtime()-start_time).gt.time_stop).or.(ii == nstep)) lstop=1
!
!	Update velocity field
!
    Call update
!
If (ifdt == 1) Call cal_dt
If (lstop == 1) goto 11
!
End Do
11 Continue
!
end_time = MPI_Wtime()
If (isavvf /= 0) Call write_restart
tot_time=end_time-start_time
!
If (nrank == 0) Write(*,*) 'Time taken: ', tot_time
!
  Call finalize_all
!
End Program cfd_incomp
