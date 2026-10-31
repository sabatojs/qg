      subroutine crank_nicol(q,dx,dy,dt,nu,N,M)

      integer N,M,i,j
      real q(N,M),dx,dy,dt,nu

      real ax,ay,a1(N),b1(N),c1(N),a2(M),b2(M),c2(M)

      real temp1(N,M),temp2(N,M)
      real rhsx(N),rhsy(M),lhsx(N),lhsy(M)

      ax=nu*dt/(dx*dx)
      ay=nu*dt/(dy*dy)

CCCCCCCCCC first direction
CCCCCCCCCCCCC j=1
      j=1
      do i=1,N
         rhsx(i) = q(i,j)
     +   +(0.*ay/4.)*(q(i,3)+q(i,1)-2.*q(i,2))
         a1(i)=-ax/2.
         b1(i)=1.+ax
         c1(i)=-ax/2.
      enddo
         call tridiag_per(a1,b1,c1,rhsx,N)
      do i=1,N
         temp1(i,j)=rhsx(i)
      enddo
CCCCCCCCCCCCC j=2,M-1
      do j=2,M-1
      do i=1,N
         rhsx(i) = q(i,j)
     +   +(ay/2.)*(q(i,j-1)+q(i,j+1)-2.*q(i,j))
         a1(i)=-ax/2.
         b1(i)=1.+ax
         c1(i)=-ax/2.
      enddo
         call tridiag_per(a1,b1,c1,rhsx,N)
      do i=1,N
         temp1(i,j)=rhsx(i)
      enddo
      enddo
CCCCCCCCCCCCC j=M
      j=M
      do i=1,N
         rhsx(i) = q(i,j)
     +   +(0.*ay/4.)*(q(i,M)+q(i,M-2)-2.*q(i,M-1))
         a1(i)=-ax/2.
         b1(i)=1.+ax
         c1(i)=-ax/2.
      enddo
         call tridiag_per(a1,b1,c1,rhsx,N)
      do i=1,N
         temp1(i,j)=rhsx(i)
      enddo

CCCCCCCCC second direction
CCCCCCCCCCCCC i=1
      i=1
      j=1
         rhsy(j) = temp1(i,j)
     +   +(ax/2.)*(temp1(2,j)+temp1(N,j)-2.*temp1(1,j))
         a2(j)=-0.*ay/4.
         b2(j)=1.+0.*ay/2.
         c2(j)=-0.*ay/4.
      do j=2,M-1
         rhsy(j) = temp1(i,j)
     +   +(ax/2.)*(temp1(2,j)+temp1(N,j)-2.*temp1(1,j))
         a2(j)=-ay/2.
         b2(j)=1.+ay
         c2(j)=-ay/2.
      enddo
      j=M
         rhsy(j) = temp1(i,j)
     +   +(ax/2.)*(temp1(2,j)+temp1(N,j)-2.*temp1(1,j))
         a2(j)=-0.*ay/4.
         b2(j)=1.+0.*ay/2.
         c2(j)=-0.*ay/4.
         call tridiag_per(a2,b2,c2,rhsy,M)
      do j=1,M
         q(i,j)=rhsy(j)
      enddo
CCCCCCCCCCCCC i=2,N-1
      do i=2,N-1
      j=1
         rhsy(j) = temp1(i,j)
     +   +(ax/2.)*(temp1(i+1,j)+temp1(i-1,j)-2.*temp1(i,j))
         a2(j)=-0.*ay/4.
         b2(j)=1.+0.*ay/2.
         c2(j)=-0.*ay/4.
      do j=2,M-1
         rhsy(j) = temp1(i,j)
     +   +(ax/2.)*(temp1(i+1,j)+temp1(i-1,j)-2.*temp1(i,j))
         a2(j)=-ay/2.
         b2(j)=1.+ay
         c2(j)=-ay/2.
      enddo
      j=M
         rhsy(j) = temp1(i,j)
     +   +(ax/2.)*(temp1(i+1,j)+temp1(i-1,j)-2.*temp1(i,j))
         a2(j)=-0.*ay/4.
         b2(j)=1.+0.*ay/2.
         c2(j)=-0.*ay/4.
         call tridiag_per(a2,b2,c2,rhsy,M)
      do j=1,M
         q(i,j)=rhsy(j)
      enddo
      enddo
CCCCCCCCCCCCC i=N
      i=N
      j=1
         rhsy(j) = temp1(i,j)
     +   +(ax/2.)*(temp1(1,j)+temp1(N-1,j)-2.*temp1(N,j))
         a2(j)=-0.*ay/4.
         b2(j)=1.+0.*ay/2.
         c2(j)=-0.*ay/4.
      do j=2,M-1
         rhsy(j) = temp1(i,j)
     +   +(ax/2.)*(temp1(1,j)+temp1(N-1,j)-2.*temp1(N,j))
         a2(j)=-ay/2.
         b2(j)=1.+ay
         c2(j)=-ay/2.
      enddo
      j=M
         rhsy(j) = temp1(i,j)
     +   +(ax/2.)*(temp1(1,j)+temp1(N-1,j)-2.*temp1(N,j))
         a2(j)=-0.*ay/4.
         b2(j)=1.+0.*ay/2.
         c2(j)=-0.*ay/4.
         call tridiag_per(a2,b2,c2,rhsy,M)
      do j=1,M
         q(i,j)=rhsy(j)
      enddo

      return
      end
