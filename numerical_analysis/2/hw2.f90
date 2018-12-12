!------------------------------------------------------
! CMSE HW 2                                      1/2017
! Computes L1, L2, and L-infinity norms
! gfortran hw2.f90 -framework Accelerate -g -fcheck=all
! jshin156@g.ucla.edu, shinjasm@msu.edu
!------------------------------------------------------


program errorAnalysis

 implicit none
 integer                        :: i, j, m
 real, dimension(5)             :: h
 real(8)                         :: PI=3.1415926535897932384626433832795, E_inf, E1, E2,E_hold,E2_hold,p
 double precision, allocatable              :: Uexact(:), U(:)

!For-whatever-h-you-want--------------

 h = (/0.1, 0.05, 0.01, 0.005, 0.001 /)
!--------------------------------------

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

! Actual Norm Calculation
! Possible future optimization: calculate U(j)-Uexact(j) store once in variable.
! Then use variable for calculation of all 3 norms.

E_inf= ABS(U(1)-Uexact(1))
 do j=2, m
    E_hold= ABS(U(j)-Uexact(j))
    if ( E_hold > E_inf) E_inf = E_hold
 enddo

do j=1, m
  E1= E1+ ABS(U(j)-Uexact(j))
enddo
E1=E1*h(i)

do j=1, m
  E2_hold= ABS(U(j)-Uexact(j))
  E2= E2+ E2_hold*E2_hold
enddo
E2=E2*h(i)
E2=E2**(0.5)

write(*,*) h(i), m, E_inf, E1, E2
call deallocateMatrix()

enddo

contains

!-------------------------------------------
!Computes the exact solution for each point , stores it in Uexact(j)
!-------------------------------------------

  subroutine initUext(N, hIN)
    implicit none
    integer, intent(IN) :: N
    real, intent(IN)    :: hIN

    allocate(Uexact(N))
   do j = 1, N
        Uexact(j) = sin(PI*j*hIN)
    enddo
  end subroutine initUext

!-------------------------------------------
! In AU=F, it calculates the F column.
! Creates the A matrix.
! Then, solves for U(estimate) using dgesv
!-------------------------------------------

  subroutine initU(N, hIN)
    implicit none
    integer, intent(IN) :: N
    real, intent(IN)    :: hIN
    double precision              :: alpha, beta
    double precision   :: a(N,N)
    integer :: ipiv(N), r
    integer :: info
    double precision, dimension(1:N) :: F

    allocate(U(N))

    alpha = 0.0
    beta = 0.0

     F(1) = -PI*PI*sin(PI*hIN) - alpha/(hIN*hIN)

     do j = 2, N-1
          F(j) = -PI*PI*sin(PI*(j)*hIN)
      !    write(*,*) F(j)  ! for testing purposes.
    enddo
    F(N) = -PI*PI*sin(PI*N*hIN) - beta/(hIN*hIN)

a(1:N,1:N)=0.0

do r =1, N
    a(r,r)= -2.0
    IF(r > 1)  a(r,r-1) = 1.0
    IF(r < N) a(r,r+1)= 1.0
enddo

!!!testing!!!
!do i= 1,N
!   do j=1,N
!       write(*,*) a(i,j)
!  enddo
!enddo

  ! does array type default to column vector?  yes
  ! dgesv computes the x in AU=F, use -framework accelerate.
  ! Answer U is stored back into F.
  call dgesv(N, 1, a, N, ipiv, F, N, info)

   do j=1, N
      U(j)=F(j)*hIN*hIN
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
