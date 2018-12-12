!------------------------------------------------------
! CMSE HW 5                                      3/2017
! Computes Steepest Descent and Conjugate Gradient Method
! gfortran hw2.f90 -framework Accelerate -g -fcheck=all
! jshin156@g.ucla.edu, shinjasm@msu.edu
!------------------------------------------------------


program SDCG

 implicit none
 integer                         :: i, j, m, ti
 real, dimension(3)              :: h
 real(8)                         :: PI=3.141592653589, E2,E2_hold, p
 double precision, allocatable   :: Uexact(:,:), U(:,:,:)
 real(8), allocatable            :: Tmatrix(:,:), Imatrix(:,:)
     real(8)    :: c

!For-whatever-h-values-you-want-----------------------
 h = (/ 0.1, 0.01/)
!-----------------------------------------------------

do i=1, size(h)

  E2=0.0
  E2_hold=0.0
  p = 1.0/h(i) -1.0

  !write(*,*) nint(p)  ! nint= nearest integer rounding otherwise int(998.9) = 998
  m=nint(p)
  !----------------------------------------------------------
  ! A^h - exists globally - overdefine my boundaries and only work within them.
  !----------------------------------------------------------
  allocate(Tmatrix(-2:m+2,-2:m+2))
  allocate(Imatrix(-2:m,-2:m+2))
Tmatrix(-2:m+2,-2:m+2) = 0.0
Imatrix(-2:m+2,-2:m+2) = 0.0

  do ti = 1, m
     Tmatrix(ti,ti)= -4.0
     Tmatrix(ti, ti-1)=1.0
     Tmatrix(ti, ti+1)=1.0
     Imatrix(ti,ti)= 1.0
  enddo




  call initUext(m,h(i))
  call initU(m,h(i))




!do j=1, m
!  E2_hold= ABS(U(j)-Uexact(j))
!  E2= E2+ E2_hold*E2_hold
!enddo
!E2=E2*h(i)
!E2=E2**(0.5)

write(*,*) h(i), m, E2
call deallocateMatrix()

enddo

contains

!-------------------------------------------------------------------------
!Computes the exact solution for each point , stores it in Uexact(j)
!-------------------------------------------------------------------------

  subroutine initUext(N, hIN)
    implicit none
    integer             :: Uext_i, ji  ! just an iterator, used here because i is used up there.
    integer, intent(IN) :: N
    real, intent(IN)    :: hIN

    allocate(Uexact(N,N))

do Uext_i = 1, N
   do ji = 1, N
        Uexact(Uext_i,j) = sin(PI*ji*hIN)*sin(PI*Uext_i*hIN)
        !write(*,*) Uexact(i,j)
   enddo
enddo

end subroutine initUext
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------


  subroutine initU(N, hIN)
    implicit none
    integer, intent(IN) :: N
    real, intent(IN)    :: hIN
    double precision, dimension(N,N,N,N)   :: A
    integer ::  i1, j1, k1, ti, tj
    real(8), dimension(N,N)   :: F
    real, dimension(1:N)      :: U1, R1
    double precision, dimension(0:N) :: alpha, beta
    double precision, dimension(0:N,N,N) :: R, W, P
    real(8)  :: tol = 1E-8

c= 0.0
allocate(U(0:N,1:N,1:N))


!----------------------------------------------------------
! U0-  Guess to be 0.
!----------------------------------------------------------
do i1 = 1, N
   do j1 = 1, N
      U(0,i1,j1) = 0.0
  enddo
enddo

!----------------------------------------------------------
! f_ij
!----------------------------------------------------------
do i1 = 1, N
    do j1 = 1, N
        F(i1,j1) = -2.0*PI*PI*sin(PI*(i1)*hIN)*sin(PI*(j1)*hIN)
          !write(*,*) F(initU_i,j)
    enddo
enddo


!----------------------------------------------------------
! R0 -- matrix multiplication on fortran?
!----------------------------------------------------------
do i1 = 1, N
   do j1 = 1, N
         call matmultiplyA(A,U,i1,0,j1)
         R(0,i1,j1) = f(i1,j1) - c
  enddo
enddo

!-------------------------------------------
! Method I - Steepest Descent
!-------------------------------------------
!k1=1
!do while ( Norm(R(k-1, 1:N, 1:N))<tol)
!  w(k1-1,1:N,1:N) = matmul(A(1:N,1:N,1:N,1:N), R(k1-1, 1:N, 1:N) )
!  alpha(k1-1) = matmul(Transpose(R(k1-1,1:N,1:N)),R(k1-1,1:N,1:N))/matmul((Transpose(R(k1-1,1:N,1:N)),w(k1-1,1:N,1:N)))
!  U(k1, 1:N, 1:N) =  U(k1-1, 1:N, 1:N) + alpha(k1-1)*R(k1-1, 1:N, 1:N)
!  R(K1, 1:N, 1:N)= R(k1-1, 1:N, 1:N) - alpha(k1-1)*w(k1-1,1:N,1:N)
!  k1=k1+1
!enddo


!-------------------------------------------
! Method II - Conjugate Gradient Method
!-------------------------------------------
!----------------------------------------------------------
! P0 -- matrix multiplication on fortran?
!----------------------------------------------------------
!do i1 = 1, N
!   do j1 = 1, N
!         P(0,i1,j1) = R(0,i1,j1) !-multiply matrix
!  enddo
!enddo


!k1=1
!do while ( Norm(R(k1, 1:N, 1:N))<tol)
!  w(k1-1,1:N,1:N) = matmul(A(1:N,1:N,1:N,1:N), P(k-1, 1:N, 1:N) )
!  alpha(k1-1) = matmul(Transpose(R(k-1,1:N,1:N)),R(k-1,1:N,1:N))/matmul((Transpose(R(k-1,1:N,1:N)),w(k1-1,1:N,1:N)))
!  U(k1, 1:N, 1:N) =  U(k1-1, 1:N, 1:N) + alpha(k1-1)*P(k1-1, 1:N, 1:N)
!  R(K, 1:N, 1:N)= R(k1-1, 1:N, 1:N) - alpha(k1-1)*w(k1-1,1:N,1:N)
! beta(k-1) =  matmul(Transpose(R(k1,1:N,1:N)),R(k1,1:N,1:N))/matmul(Transpose(R(k1-1,1:N,1:N)),R(k1-1,1:N,1:N))
!  P(K1, 1:N, 1:N)= R(K1, 1:N, 1:N)+ beta(k-1)* P(K1-1, 1:N, 1:N)
! k1=k1+1
!enddo
  end subroutine initU

!num for i1
! jout for j1

subroutine matmultiplyA(arrayA, arrayB, num, knum, jout)
   implicit none
   integer, intent(in)  :: knum, num, jout
   double precision, dimension(m,m,m,m),intent(In)   :: arrayA
   double precision, dimension(0:m,m,m), Intent(IN)   :: arrayB
   real, dimension(1:m)   :: cA ,cA1,cA2,cA3
   real, dimension(-2:m+2)    :: U31,U32,U33
   integer :: i3

   U31 = arrayB(knum,num-1,1:m)
   U32 = arrayB(knum,num,1:m)
   U33 = arrayB(knum,num+1,1:m)

cA1= matmult1(Imatrix(1:m,1:m),U31(1:m),m)
cA2= matmult1(Imatrix(1:m,1:m),U32(1:m),m)
cA3= matmult1(Imatrix(1:m,1:m),U33(1:m),m)

do i3=1,m
if(num == 1) then
     cA(i3)= cA2(i3)+cA3(i3)
else if (num .gt. 1 .or. num .lt. m) then
     cA(i3)= cA1(i3)+ cA2(i3)+cA3(i3)
else
      cA(i3)= cA1(i3)+cA2(i3)
end if

enddo
c= cA(jout)

end subroutine matmultiplyA




!-------------------------------------------
! Deallocates all the vectors
!-------------------------------------------
   subroutine deallocateMatrix()
     deallocate(Uexact)
     deallocate(U)
     deallocate(Tmatrix)
     deallocate(Imatrix)
   end subroutine deallocateMatrix
!-------------------------------------------

real function matmult1(a(:,:),b(:,1), n) result(c(:))
     integer :: i, n
     real, intent(in) :: a , b

      do i=1,n
           tmp = 0.0
          do k=1,n
              tmp = tmp + a(i,k) * b(k,1)
          enddo
           c(i) = tmp
      enddo


end function matmult1



end program SDCG







!figure out bounds of each array.
!How do I implement this for the 2D nature of my code?
