!------------------------------------------------------
! Physics 841-Project                            3/2017
! gfortran project.f90
! ./a.out > phi.dat
! Team J.
!------------------------------------------------------

program project
  implicit none
  integer                         :: i, j, k, kk, n , m, i_pole, j_pole
  real(8)                         :: dx, dy, Vmax, a, c
  real(8), allocatable            :: U(:,:,:)

  ! you can choose to divide the x and y grid with different intervals if so desired.
  dx = 15.0
  dy = 15.0

  k = 200
  Vmax = 10.0

  n = nint(1000.0/dx)
  m = nint(200.0/dy)

  i_pole= nint(500.0/dx)
  j_pole= nint(100.0/dy)

  allocate(U(0:k,0:n,0:m))

!----------------------------------------------------
! Initial guess set to zero over the whole domain
!-----------------------------------------------------

do i = 0, n
   do j = 0, m
       U(0,i,j) = 0
   enddo
enddo


do kk = 0, k

  ! set boundary for every k iteration. for the x's.
  do i = 0, N
     U(kk,i,0)= 0
     U(kk,i,m)= Vmax ! at m=20000 or 20000*dy.

  enddo

  ! set boundary for every k iteration. for the y's.
  do j = 0 , m
     U(kk,0,j)= (j*dy*Vmax)/200.0
     U(kk,N,j)= (j*dy*Vmax)/200.0
  enddo

  ! pole is grounded, so for every k iteration, at x= 100 and j up to where
  ! y = 100, the potential is set to 0.  Technically unncessary repeated below.
  do j = 0, j_pole
    U(kk, i_pole, j) = 0.0
  enddo

enddo

!---------------------------------------------
! Solver
!---------------------------------------------
a= dx*dx+dy*dy
c= (dx*dx*dy*dy*0.5)/a

do kk = 0, k-1
   do i = 1, N-1
      do j = 1, M-1
           U(kk+1,i,j) = c*((U(kk,i+1,j)+U(kk,i-1,j))/(dx*dx)+(U(kk,i,j+1)+U(kk,i,j-1))/(dy*dy))
     enddo
   enddo


! repeated pole to make double sure it's grounded every iteration.
  do j = 0, j_pole
       U(kk+1, i_pole, j) = 0.0
  enddo

enddo

do j = 0, M
    write(*, '(*(F15.5))') (U(k,i,j), i=0, N)
enddo

  deallocate(U)
end program project
