CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine fft(u,N,pwr2,itype)
      implicit none

      integer pwr2,itype,N,i
      complex u(N)
      real x(N),y(N)

      do i=1,N
         x(i)=real(u(i))
         y(i)=imag(u(i))
      enddo
      
      call sffteu( x, y, N, pwr2, itype )

      do i=1,N
         u(i)=cmplx(x(i),y(i))/(real(N)**real(itype))
      enddo

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
