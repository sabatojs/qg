      subroutine arakawa(AJ,A,B,dx,dy,N,M)
      ! calculates the arakawa jacobian partial(A,B)/partial(x,y)
      ! on an evenly spaced, doubly periodic x-y grid

      implicit none

      ! input/output
      integer N,M
      real A(N,M),B(N,M),AJ(N,M),dx,dy

      ! local workspace
      integer i,j
      real dadx(N,M),dady(N,M),dbdx(N,M),dbdy(N,M)



      ! calculate derivatives
      call d_dx(dadx,A,N,M,dx)
      call d_dx(dbdx,B,N,M,dx)
      call d_dy(dady,A,N,M,dy)
      call d_dy(dbdy,B,N,M,dy)

      ! calculate jacobian
      i=1
         j=1            
            AJ(i,j) = 0.*(
		! J1
     $		dadx(i,j)*dbdy(i,j)-dady(i,j)*dbdx(i,j)
		! J2
     $		+0.5*(a(i+1,j)*dbdy(i+1,j)-a(N,j)*dbdy(N,j))/dx 
     $		-(a(i,j+1)*dbdx(i,j+1)-a(i,j)*dbdx(i,j))/dy
		! J3
     $		-0.5*(b(i+1,j)*dady(i+1,j)-b(N,j)*dady(N,j))/dx 
     $		+(b(i,j+1)*dadx(i,j+1)-b(i,j)*dadx(i,j))/dy
     $			)/3.
         do j=2,M-1
            ! interior points
            AJ(i,j) = (
		! J1
     $		dadx(i,j)*dbdy(i,j)-dady(i,j)*dbdx(i,j)
		! J2
     $		+0.5*(a(i+1,j)*dbdy(i+1,j)-a(N,j)*dbdy(N,j))/dx 
     $		-0.5*(a(i,j+1)*dbdx(i,j+1)-a(i,j-1)*dbdx(i,j-1))/dy
		! J3
     $		-0.5*(b(i+1,j)*dady(i+1,j)-b(N,j)*dady(N,j))/dx 
     $		+0.5*(b(i,j+1)*dadx(i,j+1)-b(i,j-1)*dadx(i,j-1))/dy
     $			)/3.
	 enddo
            j=M
            AJ(i,j) = 0.*(
		! J1
     $		dadx(i,j)*dbdy(i,j)-dady(i,j)*dbdx(i,j)
		! J2
     $		+0.5*(a(i+1,j)*dbdy(i+1,j)-a(N,j)*dbdy(N,j))/dx 
     $		-(a(i,j)*dbdx(i,j)-a(i,j-1)*dbdx(i,j-1))/dy
		! J3
     $		-0.5*(b(i+1,j)*dady(i+1,j)-b(N,j)*dady(N,j))/dx 
     $		+(b(i,j)*dadx(i,j)-b(i,j-1)*dadx(i,j-1))/dy
     $			)/3.
      do i=2,N-1            
            j=1            
            AJ(i,j) = 0.*(
		! J1
     $		dadx(i,j)*dbdy(i,j)-dady(i,j)*dbdx(i,j)
		! J2
     $		+0.5*(a(i+1,j)*dbdy(i+1,j)-a(i-1,j)*dbdy(i-1,j))/dx 
     $		-(a(i,j+1)*dbdx(i,j+1)-a(i,j)*dbdx(i,j))/dy
		! J3
     $		-0.5*(b(i+1,j)*dady(i+1,j)-b(i-1,j)*dady(i-1,j))/dx 
     $		+(b(i,j+1)*dadx(i,j+1)-b(i,j)*dadx(i,j))/dy
     $			)/3.
         do j=2,M-1
            ! interior points
            AJ(i,j) = (
		! J1
     $		dadx(i,j)*dbdy(i,j)-dady(i,j)*dbdx(i,j)
		! J2
     $		+0.5*(a(i+1,j)*dbdy(i+1,j)-a(i-1,j)*dbdy(i-1,j))/dx 
     $		-0.5*(a(i,j+1)*dbdx(i,j+1)-a(i,j-1)*dbdx(i,j-1))/dy
		! J3
     $		-0.5*(b(i+1,j)*dady(i+1,j)-b(i-1,j)*dady(i-1,j))/dx 
     $		+0.5*(b(i,j+1)*dadx(i,j+1)-b(i,j-1)*dadx(i,j-1))/dy
     $			)/3.
	 enddo
            j=M
            AJ(i,j) = 0.*(
		! J1
     $		dadx(i,j)*dbdy(i,j)-dady(i,j)*dbdx(i,j)
		! J2
     $		+0.5*(a(i+1,j)*dbdy(i+1,j)-a(i-1,j)*dbdy(i-1,j))/dx 
     $		-(a(i,j)*dbdx(i,j)-a(i,j-1)*dbdx(i,j-1))/dy
		! J3
     $		-0.5*(b(i+1,j)*dady(i+1,j)-b(i-1,j)*dady(i-1,j))/dx 
     $		+(b(i,j)*dadx(i,j)-b(i,j-1)*dadx(i,j-1))/dy
     $			)/3.
      enddo
      i=N
         j=1            
            AJ(i,j) = 0.*( 
		! J1
     $		dadx(i,j)*dbdy(i,j)-dady(i,j)*dbdx(i,j)
		! J2
     $		+0.5*(a(1,j)*dbdy(1,j)-a(i-1,j)*dbdy(i-1,j))/dx 
     $		-(a(i,j+1)*dbdx(i,j+1)-a(i,j)*dbdx(i,j))/dy
		! J3
     $		-0.5*(b(1,j)*dady(1,j)-b(i-1,j)*dady(i-1,j))/dx 
     $		+(b(i,j+1)*dadx(i,j+1)-b(i,j)*dadx(i,j))/dy
     $			)/3.
         do j=2,M-1
            ! interior points
            AJ(i,j) = (
		! J1
     $		dadx(i,j)*dbdy(i,j)-dady(i,j)*dbdx(i,j)
		! J2
     $		+0.5*(a(1,j)*dbdy(1,j)-a(i-1,j)*dbdy(i-1,j))/dx 
     $		-0.5*(a(i,j+1)*dbdx(i,j+1)-a(i,j-1)*dbdx(i,j-1))/dy
		! J3
     $		-0.5*(b(1,j)*dady(1,j)-b(i-1,j)*dady(i-1,j))/dx 
     $		+0.5*(b(i,j+1)*dadx(i,j+1)-b(i,j-1)*dadx(i,j-1))/dy
     $			)/3.
	 enddo
            j=M
            AJ(i,j) = 0.*(
		! J1
     $		dadx(i,j)*dbdy(i,j)-dady(i,j)*dbdx(i,j)
		! J2
     $		+0.5*(a(1,j)*dbdy(1,j)-a(i-1,j)*dbdy(i-1,j))/dx 
     $		-(a(i,j)*dbdx(i,j)-a(i,j-1)*dbdx(i,j-1))/dy
		! J3
     $		-0.5*(b(1,j)*dady(1,j)-b(i-1,j)*dady(i-1,j))/dx 
     $		+(b(i,j)*dadx(i,j)-b(i,j-1)*dadx(i,j-1))/dy
     $			)/3.



      return
      end
