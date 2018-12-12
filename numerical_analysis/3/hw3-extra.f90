!-------------------------------------------
!    F(0) = alpha
!    F(1) =
!    do j = 2, N-1
!          F(j) = -PI*PI*sin(PI*(j)*hIN)
!      !    write(*,*) F(j)  ! for testing purposes.
!      enddo
!   F(N+1) = beta
!-------------------------------------------


!-------------------------------------------
!    F(0) = alpha
!    F(1) =
!    do j = 2, N-1
!         F(j) = -PI*PI*sin(PI*(j)*hIN)
     !    write(*,*) F(j)  ! for testing purposes.
!   enddo
!   F(N+1) = beta
!-------------------------------------------


! write(*,*) F(j)  ! for testing purposes.

set logscale x
set logscale y
f(x) = a*x+b
fit f(x) "hw2.dat" using (log($1)):(log($3)) via a, b
, "hw2.dat" using 1:4 smooth unique title 'E1', "hw2.dat" using 1:5 smooth unique title 'E2'

pause -1 "Hit any key to continue"

  !-------------------------------------------
  ! In AU=F, it calculates the F column.
  ! Creates the A matrix.
  ! Then, solves for U(estimate) using dgesv
  !-------------------------------------------

set xrange[0:1]
plot exp(x)-x+(1-2.718281828459), "method1e.dat" using 1:2 smooth unique title 'Method I'

  !-------------------------------------------
  do r =1, N
      a(r,r)= -2.0
      IF(r > 1)  a(r,r-1) = 1.0
      IF(r < N) a(r,r+1)= 1.0
  enddo
  !-------------------------------------------

Potential problems
! potential problem: j=0 index must also be filled.
! check that e^() is actually correct

! Actual Norm Calculation
! Possible future optimization: calculate U(j)-Uexact(j) store once in variable.
! Then use variable for calculation of all 3 norms.
!check bounds carefully
! potential problem: exp(x), x must be real or complex


!-------------------------------------------
!-------------------------------------------
LIST OF THINGS I STILL HAVE TO DO
!-------------------------------------------
!-------------------------------------------
CHECK DGESV!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Decide on output format
!-------------------------------------------
!-------------------------------------------
note: 0.1 is in method1-1.dat
and 0.001 is in method1-2.dat
