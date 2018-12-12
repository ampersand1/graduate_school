!-----------------------------------------------------
! Problem 1.1 from Introduction to HPC for Scientists and Engineers
!
!-Objective is to complete the fragment for a good pi approximation
!-To measure performance in MFlops/sec.
!-Lastly, estimate its latency for the operation in clock cycles.
!
!-----------------------------------------------------

PROGRAM  pi_integral

  INTEGER :: t1, t2
  DOUBLE PRECISION :: x, delta_x, sum, MFLOPS
  INTEGER, PARAMETER :: SLICES=100000000, R=10

!Vector Triad
  CALL SYSTEM_CLOCK(COUNT=t1)
  DO j=1, R
     sum = 0.d0 ; delta_x = 1.d0/SLICES
  DO i=0,SLICES-1
     x = (i+0.5)*delta_x
     sum = sum + 4.d0 / (1.d0 + x * x)
     pi = sum * delta_x
  ENDDO
  ENDDO
  CALL SYSTEM_CLOCK(COUNT=t2)


  MFLOPS = R*(SLICES-1)*2.d0/((t2-t1)*1.d6)
  WRITE(*,*) 'pi = ',pi
  WRITE(*,*)
  WRITE(*,*) 'And MFLOP/Sec is', MFLOPS
END PROGRAM pi_integral
