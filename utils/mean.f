      subroutine mean(vqyy,vq,M,P,pwr2)

CCCCC input arguments 
      integer M,P,pwr2
      real vq(M,P)

CCCCC output arguments
      real vqy(M,P),vqyy(M,P)

CCCCC local workspace
      integer j,k,M2
      real vq1(2*(M-1),P)
      real fftyr1(2*(M-1)),fftyi1(2*(M-1))
      real fftyr(4*(M-1)),fftyi(4*(M-1))
      real fftyyr(4*(M-1)),fftyyi(4*(M-1))
      complex Fy(2*(M-1),P),Fyy(2*(M-1),P)

      complex im
      real pie

      pie = 4.*atan(1.)
      im = cmplx(0.,1.)

      M2=2*(M-1)

CCCCC extend domain to twice the size (antisymmetric about y=L/2)
      do k=1,P
         vq1(1,k)=0.
         vq1(M,k)=0.
         do j=1,M-2
            vq1(j+1,k)=vq(j+1,k)
            vq1(M+j,k)=-vq(M-j,k)
         enddo
      enddo

CCCC get fourier transforms
      do k=1,P
         do j=1,M2
            fftyr1(j)=vq1(j,k)
            fftyi1(j)=0.
         enddo
         call sffteu(fftyr1,fftyi1,M2,pwr2,1)
         do j=1,M2
            ll=(j-1.)*2.*pie/L
            if (j.gt.M2/2) ll=(-M2+j-1.)*2.*pie/L
            Fy(j,k)=im*ll*cmplx(fftyr1(j),fftyi1(j))
            Fyy(j,k)=-ll*ll*cmplx(fftyr1(j),fftyi1(j))
         enddo
      enddo


CCCCC transform to physical space
      do k=1,P
         do j=1,2*M2
            fftyr(j)=real(Fy(j,k))
            fftyi(j)=0.*imag(Fy(j,k))
            fftyyr(j)=0.*real(Fyy(j,k))
            fftyyi(j)=imag(Fyy(j,k))
         enddo
         call sffteu(fftyr,fftyi,2*M2,pwr2*2,-1)
         call sffteu(fftyyr,fftyyi,2*M2,pwr2*2,-1)
         do j=1,M2,2
            vqy((j+1)/2,k) = 2.*fftyr(j)
            vqyy((j+1)/2,k) = 2.*fftyyr(j)
         enddo
      enddo

      return
      end
