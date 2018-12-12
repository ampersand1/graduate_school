!------------------------------------------------------
! CMSE HW 3                                      2/2017
! Computes L1, L2, and L-infinity norms
! look at lecture 4 and lecture 7 notes
! gfortran hw2.f90 -framework Accelerate -g -fcheck=all
! jshin156@g.ucla.edu, shinjasm@msu.edu
!------------------------------------------------------


program errorAnalysis

 implicit none
 integer                         :: i, j, m
 real, dimension(5)              :: h
 real(8)                         :: PI=3.141592653589, E1, E2, E_inf,E_hold,E2_hold, p
 real(8)                         :: e= 2.718281828459
 double precision, allocatable   :: Uexact(:), U(:)

!For-whatever-h-values-you-want-----------------------
 h = (/ 0.1, 0.05, 0.01, 0.005, 0.001 /)
!-----------------------------------------------------

do i=1, size(h)

  E1=0.0
  E2=0.0
  E_hold=0.0
  E2_hold=0.0
  E_inf=0.0
  p = 1.0/h(i) -1.0

  !write(*,*) nint(p)  ! nint= nearest integer rounding otherwise int(998.9) = 998
  m=nint(p)
  call initUext(m,h(i))
  call initU(m,h(i))


E_inf= ABS(U(1)-Uexact(1))
 do j=2, m
    E_hold= ABS(U(j)-Uexact(j))
    if ( E_hold > E_inf) E_inf = E_hold
 enddo

do j=1, m
  E2_hold= ABS(U(j)-Uexact(j))
  E2= E2+ E2_hold*E2_hold
enddo
E2=E2*h(i)
E2=E2**(0.5)

write(*,*) h(i), m, E_inf, E2
call deallocateMatrix()

enddo

contains

!-------------------------------------------------------------------------
!Computes the exact solution for each point , stores it in Uexact(j)
!-------------------------------------------------------------------------

  subroutine initUext(N, hIN)
    implicit none
    integer, intent(IN) :: N
    real, intent(IN)    :: hIN

    allocate(Uexact(0:N+1))
   do j = 0, N+1
        Uexact(j) = exp(hIN*j)-j*hIN+(1-e)
  enddo
  end subroutine initUext

!-------------------------------------------------------------------------


  subroutine initU(N, hIN)
    implicit none
    integer, intent(IN) :: N
    real, intent(IN)    :: hIN
    double precision              :: alpha, beta
    double precision   :: a(0:N+1,0:N+1)
    integer :: ipiv(N+2), r
    integer :: info
    double precision, dimension(0:N+1) :: F

    allocate(U(0:N+1))

!  These changes based on the problem.
    alpha = 0.0  ! sigma according to notes, We're talking Neuman BP
    beta = 0.0

!-------------------------------------------
! Method I
!-------------------------------------------
!  F(0) = alpha
!    do j = 1, N
!          F(j) = exp(hIN*j)
!    enddo
!    F(N+1) = beta

!a(0:N+1,0:N+1)=0.0

!a(0,0)= -hIN
!a(0,1)=  hIN
!a(N+1,N+1)= hIN*hIN
!do r =1, N
!    a(r,r)= -2.0
!    a(r,r-1) = 1.0
!    a(r,r+1)= 1.0
!enddo

!-------------------------------------------
! Method II
!-------------------------------------------
!   F(0) = alpha + (hIN/2.0)*exp(0.0)
!     do j = 1, N
!          F(j) = exp(hIN*j)
!    enddo
!    F(N+1) = beta

!a(0:N+1,0:N+1)=0.0

!a(0,0)= -hIN
!a(0,1)=  hIN
!a(N+1,N+1)= hIN*hIN
!do r =1, N
!    a(r,r)= -2.0
!    a(r,r-1) = 1.0
!    a(r,r+1)= 1.0
!enddo

!-------------------------------------------
! Method III
!-------------------------------------------
F(0) = alpha
     do j = 1, N
          F(j) = exp(hIN*j)
     enddo
  F(N+1) = beta

a(0:N+1,0:N+1)=0.0

a(0,0)= 1.5*hIN
a(0,1)= -2.0*hIN
a(0,2)= 0.5*hIN
a(N+1,N+1)= hIN*hIN
do r =1, N
    a(r,r)= -2.0
    a(r,r-1) = 1.0
    a(r,r+1)= 1.0
enddo

!-------------------------------------------

!This might screw up
  call dgesv(N+2, 1, a, N+2, ipiv, F, N+2, info)

   do j=0, N+1
      U(j)=F(j)*hIN*hIN
!     write(*,*) j*hIN, U(j)
   enddo

  end subroutine initU

!-------------------------------------------
! Deallocates all the vectors
!-------------------------------------------
   subroutine deallocateMatrix()
     deallocate(Uexact)
     deallocate(U)
   end subroutine deallocateMatrix
!-------------------------------------------




end program errorAnalysis
