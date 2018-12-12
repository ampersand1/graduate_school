!------------------------------------------
! Finite Difference Approximation
! U(x) = x^3/2 and X-bar=0
! jshin156@g.ucla.edu, shinjasm@msu.edu
!-------------------------------------------


program fda
 implicit none
 integer                        :: i
 double precision               :: rfda, u, x, u_F, uh, uh_F, dplus, &
                                   uprime, uprime_F
 real, dimension(5) :: h

!For-whatever-X-you-want--------------
 x = 0.0
 h = (/ 1e-1, 0.01, 0.001, 0.0001, 0.00001 /)
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


! all the derivative components here

function u_F(x) result(r)
  double precision, intent(in) :: x
  double precision :: r
  r= x**1.5
end function u_F


function uh_F(x,h) result(r)
  double precision, intent(in) :: x
  real, intent(in) :: h
  double precision :: r
  r= (x+h)**1.5
end function uh_F

function uprime_F(x) result(r)
  double precision, intent(in) :: x
  double precision :: r
  r= 1.5*x**0.5
end function uprime_F
