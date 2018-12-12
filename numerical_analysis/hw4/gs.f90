!------------------------------------------------------
! CMSE HW 4                                      2/2017
! I) Jacobi Method
! II) G-S Method
! III) SOR
! gfortran hw4.f90 -framework Accelerate -g -fcheck=all
! jshin156@g.ucla.edu, shinjasm@msu.edu
!------------------------------------------------------

program errorAnalysis
  implicit none
  integer                         :: i, j, k, m, niter, ii,jj, k_loop
  real, dimension(1)              :: h
  integer, dimension(5)           :: karray
  real(8)                         :: PI=3.141592653589, E2, E2_hold, p
  real(8), allocatable            :: Uexact(:,:), U(:,:,:)


! Set your h's to whatever you desire-----------------------
  h = (/0.1/)
  karray = (/ 1, 3, 5, 7, 9/)

  !h=(/0.01/)
  !karray = (/ 1, 50, 100, 150, 200/)

  !h=(/0.001/)
  !karray = (/ 1, 100, 500, 1000, 1500/)
! ----------------------------------------------------------

do i=1, size(h)
  do k_loop=1,size(karray)

  E2=0.0
  E2_hold=0.0

  p = 1.0/h(i) -1.0
  m=nint(p)


  niter = karray(k_loop)

  call initUext(m,h(i))
  call initU(m,h(i),niter)


  do ii=1, m
    do jj=1, m
     E2_hold = ABS(U(niter,ii,jj)-Uexact(ii,jj))
     E2 = E2+ E2_hold*E2_hold
    enddo
  enddo
  E2=E2*h(i)
  E2=E2**(0.5)


  write(*,*) h(i), m , karray(k_loop), E2
  call deallocateMatrix()

enddo
enddo

contains

!-------------------------------------------------------------------------
!Computes the exact solution for each point , stores it in Uexact(j)
!-------------------------------------------------------------------------

  subroutine initUext(N, hIN)
    implicit none
    integer             :: Uext_i  ! just an iterator, used here because i is used up there.
    integer, intent(IN) :: N
    real, intent(IN)    :: hIN

    allocate(Uexact(N,N))

do Uext_i = 1, N
   do j = 1, N
        Uexact(Uext_i,j) = sin(PI*j*hIN)*sin(PI*Uext_i*hIN)
        !write(*,*) Uexact(i,j)
   enddo
enddo

end subroutine initUext
!-------------------------------------------------------------------------


  subroutine initU(N, hIN, n_it)
    implicit none
    integer             :: initU_i
    integer, intent(IN) :: N, n_it
    real, intent(IN)    :: hIN
    double precision, dimension(N,N) :: F

    allocate(U(0:n_it,0:N+1,0:N+1))

!----------------------------------------------------
! Initial guess set to zero over the whole domain
!-----------------------------------------------------

do initU_i = 0, N+1
  do j = 0, N+1
      U(0,initU_i,j) = 0
      !write(*,*) U(0,initU_i,j)
  enddo
enddo

! Set the boundary to 0 for all K.

do k = 0, n_it
  do initU_i = 0, N+1
   U(k,initU_i,0)= 0
   U(k,initU_i,N+1)= 0
  enddo

  do j = 0 , N+1
   U(k,0,j)= 0
   U(k,N+1,j)= 0
 enddo
enddo

!----------------------------------------------------------
! f_ij
!----------------------------------------------------------
do initU_i = 1, N
   do j = 1, N
      F(initU_i,j) = -2.0*PI*PI*sin(PI*(initU_i)*hIN)*sin(PI*(j)*hIN)
      !write(*,*) F(initU_i,j)
   enddo
enddo

!-------------------------------------------
! Method II - G-S  !!!!! Do I need to set extra initial conditions?
!-------------------------------------------
do k = 0, n_it-1
   do initU_i = 1, N
      do j = 1, N
         U(k+1,initU_i,j) = 0.25*(u(k+1,initU_i-1,j)+u(k+1,initU_i,j-1)+u(k,initU_i+1,j)+u(k,initU_i,j+1))-0.25*hIN*hIN*F(initU_i,j)
      enddo
  enddo
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
