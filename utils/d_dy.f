CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine d_dy(dAdy,A,N,M,dy)
                                                                                                                                                 
      integer N,M
      real A(N,M),dAdy(N,M)

      integer i,j                                                                                                                                                 
                                                                                                                                                 
      do i=1,N
         dAdy(i,1) = (A(i,2)-A(i,1))/dy
      do j=2,M-1
         dAdy(i,j) = 0.5*(A(i,j+1)-A(i,j-1))/dy
      enddo             
         dAdy(i,M) = (A(i,M)-A(i,M-1))/dy
      enddo
                                                                                                                                    
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

