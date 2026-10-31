CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine d_dx(dAdx,A,N,M,dx)
                                                                                                                                                 
      integer N,M
      real A(N,M),dAdx(N,M)

      integer i,j                                                                                                                                                 
                                                                                                                                                 
      do j=1,M
         dAdx(1,j) = 0.5*(A(2,j)-A(N,j))/dx
      do i=2,N-1
         dAdx(i,j) = 0.5*(A(i+1,j)-A(i-1,j))/dx
      enddo             
         dAdx(N,j) = 0.5*(A(1,j)-A(N-1,j))/dx
      enddo
                                                                                                                                    
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

