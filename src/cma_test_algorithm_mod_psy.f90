!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
program invoke_12_columnwise_op_mul_kernel_type

use constants_mod, only: r_def, r_double
use columnwise_op_mul_kernel_mod, only: columnwise_op_mul_kernel_code
use dino_mod, only: dino_type
!$ use omp_lib

implicit none 

integer :: cell
integer :: last_halo_cell
integer :: ncell_2d

integer :: mm_vel_v_inv_diag_cma_nrow, mm_vel_v_inv_diag_cma_ncol,           &
           mm_vel_v_inv_diag_cma_bandwidth, mm_vel_v_inv_diag_cma_alpha,     &
           mm_vel_v_inv_diag_cma_beta, mm_vel_v_inv_diag_cma_gamma_m,        &
           mm_vel_v_inv_diag_cma_gamma_p
real(kind=r_def), allocatable :: mm_vel_v_inv_diag_cma_matrix(:,:,:)

integer :: grad_v_cma_nrow, grad_v_cma_ncol, grad_v_cma_bandwidth,           &
           grad_v_cma_alpha, grad_v_cma_beta, grad_v_cma_gamma_m,            &
           grad_v_cma_gamma_p

real(kind=r_def), allocatable :: grad_v_cma_matrix(:,:,:)

integer :: mm_vel_v_inv_grad_v_cma_nrow, mm_vel_v_inv_grad_v_cma_ncol,       &
           mm_vel_v_inv_grad_v_cma_bandwidth, mm_vel_v_inv_grad_v_cma_alpha, &
           mm_vel_v_inv_grad_v_cma_beta, mm_vel_v_inv_grad_v_cma_gamma_m,    &
           mm_vel_v_inv_grad_v_cma_gamma_p
real(kind=r_def), allocatable :: mm_vel_v_inv_grad_v_cma_matrix(:,:,:)

real(kind=r_def), allocatable :: known_output(:,:,:)

integer, parameter :: attempts = 100
!$ real(kind=r_double) :: tinit, tend
!$ real(kind=r_double) :: timings(attempts)
type(dino_type) :: dino_in
integer :: i, j, k, iter
real(kind=r_def) :: delta
real(kind=r_def), parameter :: tolerance = 1e-12
logical :: comparison_ok = .true.
logical, parameter :: debug = .false.

dino_in = dino_type()

write(6,'(A)') 'Reading dinodump'

call dino_in%input_scalar(last_halo_cell)
call dino_in%input_scalar(ncell_2d)
call dino_in%input_scalar(mm_vel_v_inv_diag_cma_nrow)
call dino_in%input_scalar(mm_vel_v_inv_diag_cma_ncol)
call dino_in%input_scalar(mm_vel_v_inv_diag_cma_bandwidth)
call dino_in%input_scalar(mm_vel_v_inv_diag_cma_alpha)
call dino_in%input_scalar(mm_vel_v_inv_diag_cma_beta)
call dino_in%input_scalar(mm_vel_v_inv_diag_cma_gamma_m)
call dino_in%input_scalar(mm_vel_v_inv_diag_cma_gamma_p)
call dino_in%input_scalar(grad_v_cma_nrow)
call dino_in%input_scalar(grad_v_cma_ncol)
call dino_in%input_scalar(grad_v_cma_bandwidth)
call dino_in%input_scalar(grad_v_cma_alpha)
call dino_in%input_scalar(grad_v_cma_beta)
call dino_in%input_scalar(grad_v_cma_gamma_m)
call dino_in%input_scalar(grad_v_cma_gamma_p)
call dino_in%input_scalar(mm_vel_v_inv_grad_v_cma_nrow)
call dino_in%input_scalar(mm_vel_v_inv_grad_v_cma_ncol)
call dino_in%input_scalar(mm_vel_v_inv_grad_v_cma_bandwidth)
call dino_in%input_scalar(mm_vel_v_inv_grad_v_cma_alpha)
call dino_in%input_scalar(mm_vel_v_inv_grad_v_cma_beta)
call dino_in%input_scalar(mm_vel_v_inv_grad_v_cma_gamma_m)
call dino_in%input_scalar(mm_vel_v_inv_grad_v_cma_gamma_p)

if (debug) then
  write(6, *) last_halo_cell
  write(6, *) ncell_2d
  write(6, *) mm_vel_v_inv_diag_cma_nrow
  write(6, *) mm_vel_v_inv_diag_cma_ncol
  write(6, *) mm_vel_v_inv_diag_cma_bandwidth
  write(6, *) mm_vel_v_inv_diag_cma_alpha
  write(6, *) mm_vel_v_inv_diag_cma_beta
  write(6, *) mm_vel_v_inv_diag_cma_gamma_m
  write(6, *) mm_vel_v_inv_diag_cma_gamma_p
  write(6, *) grad_v_cma_nrow
  write(6, *) grad_v_cma_ncol
  write(6, *) grad_v_cma_bandwidth
  write(6, *) grad_v_cma_alpha
  write(6, *) grad_v_cma_beta
  write(6, *) grad_v_cma_gamma_m
  write(6, *) grad_v_cma_gamma_p
  write(6, *) mm_vel_v_inv_grad_v_cma_nrow
  write(6, *) mm_vel_v_inv_grad_v_cma_ncol
  write(6, *) mm_vel_v_inv_grad_v_cma_bandwidth
  write(6, *) mm_vel_v_inv_grad_v_cma_alpha
  write(6, *) mm_vel_v_inv_grad_v_cma_beta
  write(6, *) mm_vel_v_inv_grad_v_cma_gamma_m
  write(6, *) mm_vel_v_inv_grad_v_cma_gamma_p
end if
! Allocate Matrices
allocate(mm_vel_v_inv_diag_cma_matrix(mm_vel_v_inv_diag_cma_bandwidth,&
                                      mm_vel_v_inv_diag_cma_nrow,&
                                      ncell_2d))

allocate(grad_v_cma_matrix(grad_v_cma_bandwidth, &
                           grad_v_cma_nrow,&
                           ncell_2d))

allocate(known_output(mm_vel_v_inv_grad_v_cma_bandwidth,&
                      mm_vel_v_inv_grad_v_cma_nrow,&
                      ncell_2d))

! populate matrices
!

call dino_in%input_array(mm_vel_v_inv_diag_cma_matrix, &
                         mm_vel_v_inv_diag_cma_bandwidth, &
                         mm_vel_v_inv_diag_cma_nrow,&
                         ncell_2d)

call dino_in%input_array(grad_v_cma_matrix,&
                         grad_v_cma_bandwidth, &
                         grad_v_cma_nrow,&
                         ncell_2d)

call dino_in%input_array(known_output,&
                         mm_vel_v_inv_grad_v_cma_bandwidth,&
                         mm_vel_v_inv_grad_v_cma_nrow,&
                         ncell_2d)

call dino_in%io_close()
write(6,'(A)') 'Dinodump read complete'

if (debug) then
  write(6,*) 'Max, min A = ', maxval(mm_vel_v_inv_diag_cma_matrix), minval(mm_vel_v_inv_diag_cma_matrix)
  write(6,*) 'Max, min B = ', maxval(grad_v_cma_matrix), minval(grad_v_cma_matrix)
  write(6,*) 'MAX, min, KGO', maxval(known_output), minval(known_output)
end if


allocate(mm_vel_v_inv_grad_v_cma_matrix(mm_vel_v_inv_grad_v_cma_bandwidth,&
                                        mm_vel_v_inv_grad_v_cma_nrow,&
                                        ncell_2d))

!
! Call kernel
!
write(6,'(A,I0,A,I0,A)') 'Beginning computation over ', last_halo_cell, ' cells ', attempts, ' times' 
do iter=1,attempts
  mm_vel_v_inv_grad_v_cma_matrix(:,:,:) = 0.0
!$ tinit = omp_get_wtime()
  do cell=1,last_halo_cell
    !
    call columnwise_op_mul_kernel_code(cell, ncell_2d, mm_vel_v_inv_diag_cma_matrix, mm_vel_v_inv_diag_cma_nrow, &
       mm_vel_v_inv_diag_cma_ncol, mm_vel_v_inv_diag_cma_bandwidth, mm_vel_v_inv_diag_cma_alpha, mm_vel_v_inv_diag_cma_beta, &
       mm_vel_v_inv_diag_cma_gamma_m, mm_vel_v_inv_diag_cma_gamma_p, grad_v_cma_matrix, grad_v_cma_nrow, grad_v_cma_ncol, &
       grad_v_cma_bandwidth, grad_v_cma_alpha, grad_v_cma_beta, grad_v_cma_gamma_m, grad_v_cma_gamma_p, mm_vel_v_inv_grad_v_cma_matrix, &
       mm_vel_v_inv_grad_v_cma_nrow, mm_vel_v_inv_grad_v_cma_ncol, mm_vel_v_inv_grad_v_cma_bandwidth, mm_vel_v_inv_grad_v_cma_alpha, &
       mm_vel_v_inv_grad_v_cma_beta, mm_vel_v_inv_grad_v_cma_gamma_m, mm_vel_v_inv_grad_v_cma_gamma_p)
  end do
!$ tend = omp_get_wtime()
!$ timings(iter) = tend-tinit
end do
!$ write(6,'(A, F16.4)') 'Mean Time taken = ', SUM(timings) / REAL(attempts)
write(6,'(A)') 'End of computation'
if (debug) then
  write(6,*) 'MAX, min, output', maxval(mm_vel_v_inv_grad_v_cma_matrix), minval(mm_vel_v_inv_grad_v_cma_matrix)
end if
! check the answer
do k = 1, ncell_2d
  do j = 1, mm_vel_v_inv_grad_v_cma_nrow
    do i = 1, mm_vel_v_inv_grad_v_cma_bandwidth
      delta = abs(mm_vel_v_inv_grad_v_cma_matrix(i,j,k) - known_output(i,j,k))
      if(delta > tolerance) then
        comparison_ok = .false.
        write(6, *) i,j,k,mm_vel_v_inv_grad_v_cma_matrix(i,j,k), known_output(i,j,k), delta
      end if
    end do
  end do
end do

if (comparison_ok) then
  write(6,'(A)') 'Output from kernel compares'
else
  write(6,'(A)') 'FAIL: Output from kernel does not compare'
end if

deallocate(mm_vel_v_inv_grad_v_cma_matrix)
deallocate(grad_v_cma_matrix)
deallocate(known_output)
deallocate(mm_vel_v_inv_diag_cma_matrix)

end program invoke_12_columnwise_op_mul_kernel_type
