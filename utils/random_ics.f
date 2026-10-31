      subroutine random_ics(q,N,M,P)

      implicit none

C input
      integer N,M,P

C output
      real q(N,M,P)

C local variables
      integer randseed
      integer i,j,k
      real qs(N,M-1)
      real pie

C external functions
      real ran1   


      pie = 4.*atan(1.)
      q = 0.

      randseed = -23
      call initial(N,M-1,qs,randseed)
      do i=1,N
            q(i,1,1)=0.
         do j=2,M-1
            q(i,j,1)=qs(i,j)*sin((j-1.)*pie/M)*sin((j-1.)*pie/M)
Cexp(-(j-M/2+0.5)*(j-M/2+0.5)*10./(M*M))
         enddo
            q(i,M,1)=0.
      enddo
      randseed = -47
      call initial(N,M-1,qs,randseed)
      do i=1,N
            q(i,1,P)=0.
         do j=2,M-1
            q(i,j,P)=qs(i,j)*sin((j-1.)*pie/M)*sin((j-1.)*pie/M)
Cexp(-(j-M/2+0.5)*(j-M/2+0.5)*10./(M*M))
         enddo
            q(i,M,P)=0.
      enddo

      return
      end
