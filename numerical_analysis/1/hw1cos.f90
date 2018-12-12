!------------------------------------------
! Finite Difference Approximation
! U(x) = x^3/2 and X-bar=0
!
!-------------------------------------------


program fda
 implicit none
 integer                        :: i
 double precision               :: rfda, u, u_F, uh, uh_F, dplus, &
                                   uprime, uprime_F
 real, dimension(5) :: h
 real :: x

!For-whatever-X-you-want--------------
 x = 1.0
 h = (/ 0.1, 0.01, 0.001, 0.0001, 0.00001 /)
!--------------------------------------

do i=1, size(h)
 u =  u_F(x)
 uh = uh_F(x,h(i))
 dplus=(uh-u)/h(i)
 uprime=uprime_F(x)
 rfda= dplus-uprime
 write(*,*) h(i), rfda
enddo

end program fda


function u_F(x) result(r)
  real, intent(in) :: x
  double precision :: r
  r= cos(x)
end function u_F


function uh_F(x,h) result(r)
  real, intent(in) :: x, h
  double precision :: r
  r= cos(x+h)
end function uh_F

function uprime_F(x) result(r)
  real, intent(in) :: x
  double precision:: r
  r= -sin(x)
end function uprime_F
