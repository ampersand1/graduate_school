!------------------------------------------
! HW 2
!
! jshin156@g.ucla.edu, shinjasm@msu.edu
!-------------------------------------------


program errorAnalysis
 implicit none
 integer                        :: i, j, m
 real, dimension(5)             :: h
 real                           :: PI=3.14159265, E_inf, E1, E2,E_hold
 real, allocatable              :: Uexact(:), U(:), F(:)


!For-whatever-X-you-want--------------

 h = (/ 0.1, 0.05, 0.01, 0.005, 0.001 /)
!--------------------------------------

do i=1, size(h)

  E1=0.0
  E2=0.0

   m = 1/h(i) -1
  call initUext(m,h(i))
  call initU(m,h(i))

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
  E2= E2+ ABS(U(j)-Uexact(j))
  E2=E2*E2
enddo
E2=E2*h(i)
E2=E2**(0.5)

!write(*,*) h(i), m, E_inf, E1, E2
call deallocateMatrix()

enddo


contains
!-------------------------------------------
  subroutine initUext(N, hIN)
    implicit none
    integer, intent(IN) :: N
    real, intent(IN)    :: hIN

    allocate(Uexact(N))
   do j = 1, N
        Uexact(j) = sin(PI*j*hIN)
        !write(*,*) Uexact(j)
    enddo
  end subroutine initUext

!-------------------------------------------

  subroutine initU(N, hIN)
    implicit none
    integer, intent(IN) :: N
    real, intent(IN)    :: hIN
    real                :: alpha, beta

    allocate(U(N))
    allocate(F(N))

    alpha = 0.0
    beta = 0.0

     F(1) = -PI*PI*sin(PI*hIN) - alpha/(hIN*hIN)

     do j = 2, N-1
          F(j) = -PI*PI*sin(PI*(j-1)*hIN)
          write(*,*) F(j)  ! for testing purposes.
    enddo

    F(N) = -PI*PI*sin(PI*N*hIN) - beta/(hIN*hIN)


  end subroutine initU


   subroutine deallocateMatrix()
     deallocate(Uexact)
     deallocate(U)
     deallocate(F)
   end subroutine deallocateMatrix
!-------------------------------------------

end program errorAnalysis


!-------------------------------------------
    !This is for U(0), but for fortran indexing purposes, called it "alpha"
    !Defined this way to make it easier to change function if so desired.
  !  alpha = 0.0
  !  U(1) = 0.0
  !  U(2) = -PI*PI*sin(PI*hIN)*hIN*hIN - alpha + 2*U(1)

  ! do j = 3, N
  !      holdingNo = -PI*PI*sin(PI*(j-1)*hIN)
  !      U(j) = hIN*hIN*holdingNo - U(j-1)+ 2*U(j)
  !      write(*,*) U(j)  ! for testing purposes.
  !  enddo
!-------------------------------------------
