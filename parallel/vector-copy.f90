!-----------------------------------------------
!vector-Copy.f90
!
!This code copies over one vector into another
!and then returns the memory bandwidth of the code in MB/sec.
!
!-jshin
!-----------------------------------------------


program vectorCopy
    implicit none
    integer :: i, k
    integer(8), parameter :: n = 2000000
    integer, parameter :: nk = 10000    ! nk is the number of trials
    double precision, dimension(0:n-1) :: A, C
    double precision :: r,E,S,bandwidth


  !This loop fills in a vector.
  do i=0, n, 1
        ! a random number between 0 and 100
        call random_number(r)
        a(i) = r*100
  enddo

!This loop does the actual copying.
!Also, it calculates the time it took before and after.
call get_walltime(S)

  Do k=0, nk , 1
      Do i=0, n, 1
          c(i)=a(i)
      enddo
      !The dummy function is to trick compiler optimizations.
      call dummy(c(i))
  enddo

call get_walltime(E)


 bandwidth = 16*nk*n/((E-S)*1.d6)
 write(*,*) 'Bandwidth is', bandwidth


end program vectorCopy
