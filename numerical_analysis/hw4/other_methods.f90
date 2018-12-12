
!-------------------------------------------
! Method II - G-S  !!!!! Do I need to set extra initial conditions?
!-------------------------------------------
!do k = 0, niter
!   do i = 1, N
!      do j = 1, N
!         U(k+1,i,j) = 0.25*(u(k+1,i-1,j)+u(k+1,i,j-1)+u(k,i+1,j)+u(k,i,j+1))-0.5*hIN*hIN*F(i,j)
!      enddo
!  enddo
!enddo

!-------------------------------------------
! Method III - SOR
!-------------------------------------------
!do k = 0, niter
!   do i = 1, N
!      do j = 1, N
!         U(k+1,i,j) = u(k,i,j)+ w*(u(k+1,i-1,j)+u(k+1,i,j-1)+u(k,i+1,j)+u(k,i,j+1))-0.5*hIN*hIN*F(i,j)-u(k,i,j))
!      enddo
!  enddo
!enddo
!-------------------------------------------
