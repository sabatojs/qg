      subroutine four1(data,nn,isign)
      integer isign,nn
      real data(2*nn)
C
C
      integer i, istep,j,m,mmax,n
      real tempi,tempr
      double precision theta,wi,wpi,wpr,wr,wtemp

      n=2*nn
      j=1
      do i=1,n,2
         if (j.gt.i) then
            tempr=data(j)
            tempi=data(j+1)
            data(j)=data(i)
            data(j+1)=data(i+1)
            data(i)=tempr
            data(i+1)=tempi
         endif
         m=nn
  1      if ((m.ge.2).and.(j.gt.m)) then
            j=j-m
            m=m/2
         goto 1
         endif
         j=j+m
      enddo
      mmax=2
  2   if (n.gt.mmax) then
         istep=2*mmax
         theta=6.28318530717959d0/(isign*mmax)
         wpr=-2.d0*sin(0.5d0*theta)**2
         wpi=sin(theta)
         wr=1.d0
         wi=0.d0
         do m=1,mmax,2
            do i=m,n,istep
               j=i+mmax
               tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
               tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
               data(j)=data(i)-tempr
               data(j+1)=data(i+1)-tempi
               data(i)=data(i)+tempr
               data(i+1)=data(i+1)+tempi
            enddo
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
         enddo
         mmax=istep
      goto 2
      endif
      return
      end
