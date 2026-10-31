      subroutine stream(psi,v,u,q,qy,
     $			dy,dz,eps,Lx,Ly,N,M,P,pwr2,nlin,k0)

CCCC input arguments
      integer N,M,P,pwr2,k0
      real q(N,M,P),qy(M,P)
      real dy,dz,eps,Lx,Ly,nlin

CCCC output arguments
      real psi(N,M,P),v(N,M,P),u(M,P)

CCCC local workspace
      integer i,j,k
      real kk,ll
      real q1(N,M-1,P),qy1(M-1,P)
      complex F(N,M-1,P)
      real Fy(M-1,P)
      complex phi(2*N,M-1,P),vsp(2*N,M-1,P),usp(M-1,P)
      real tmp1(N,M-1),tmp2(N,M-1),tmp3(N,M-1),tmp4(N,M-1)
      complex a(P-1),b(P),c(P-1),rhs(P)
      real fftxr1(N),fftxi1(N)
      real fftxr2(N),fftxi2(N)
      real fftyr1(M-1),fftyi1(M-1)
      real fftyr2(M-1),fftyi2(M-1)
      real fftxr(2*N),fftxi(2*N),fftxvr(2*N),fftxvi(2*N)
      real fftyr(2*(M-1)),fftyi(2*(M-1))
      real fftyrp(2*(M-1)),fftyip(2*(M-1))
      real fftyrv(2*(M-1)),fftyiv(2*(M-1))
      real ubtsp(M-1),phibt(2*N,M-1)
      real Fybt(M-1),Fbt(N,M-1)
      real qybt(M-1),qbt(N,M-1)
      complex im
      real pie
      integer k1,k2

      pie = 4.*atan(1.)
      im = cmplx(0.,1.)

      phisp = 0.
      vsp = 0.
      ubtsp = 0.
      usp = 0.

      psi=0.
      v=0.
      u=0.
      

CCCCC if it's linear just do that wavenumber
      if (nlin.eq.1.) then
         k1=2
         k2=N/2
      else
         k1=k0+1
         k2=k0+1
      endif

      do j=1,M-1
            qybt(j)=0.
         do k=1,P
            qybt(j)=((k-1.)*qybt(j) + qy(j,k))/k
         enddo
C            do i=1,N
C                  qbt(i,j)=0.
C               do k=1,P
C                  qbt(i,j)=((k-1.)*qbt(i,j) + q(i,j,k))/k
C               enddo         
C            enddo
      enddo

CCCC get fourier transforms
CCCC sin transform in y for mean barotropic flow
         do j=1,M-1
            fftyr1(j)=-qybt(j)
         enddo
         call sinft(fftyr1,M-1)                
         do j=1,M-1
            Fybt(j)=fftyr1(j)
         enddo

CCCC sin transform in y for barotropic eddies
C         do i=1,N
C            do j=1,M-1
C               fftyr1(j)=qbt(i,j)
C            enddo
C            call sinft(fftyr1,M-1)                
C            do j=1,M-1
C               tmp4(i,j)=fftyr1(j)
C            enddo
C         enddo
CCCC full fft in x for barotropic eddies
C         do j=1,M-1
C            do i=1,N
C               fftxr1(i)=tmp4(i,j)
C               fftxi1(i)=0.
C            enddo
C               call sffteu(fftxr1,fftxi1,N,pwr2,1)       
C            do i=1,N
C               Fbt(i,j)=cmplx(fftxr1(i),fftxi1(i))
C            enddo
C         enddo

      do k=1,P
CCCC sin transform in y for mean
         do j=1,M-1
            fftyr1(j)=-qy(j,k) + qybt(j)
         enddo
         call sinft(fftyr1,M-1)                
         do j=1,M-1
            Fy(j,k)=fftyr1(j)
         enddo

CCCC sin transform in y for eddies
         do i=1,N
            do j=1,M-1
               fftyr1(j)=q(i,j,k) 
C               fftyr1(j)=q(i,j,k) - qbt(i,j)
            enddo
            call sinft(fftyr1,M-1)                
            do j=1,M-1
               tmp4(i,j)=fftyr1(j)
            enddo
         enddo
CCCC full fft in x for eddies
         do j=1,M-1
            do i=1,N
               fftxr1(i)=tmp4(i,j)
               fftxi1(i)=0.
            enddo
               call sffteu(fftxr1,fftxi1,N,pwr2,1)       
            do i=1,N
               F(i,j,k)=cmplx(fftxr1(i),fftxi1(i))
            enddo
         enddo
      enddo

CCCCC do eddies
      do i=k1,k2
         kk=(i-1.)*2.*pie/Lx
         do j=1,M-1
               ll=(j-1.)*pie/Ly
               b(1)=-(ll*ll + kk*kk +(eps/(dz*dz)))
               c(1)=eps/(dz*dz)
               rhs(1)=F(i,j,1)
            do k=2,P-1
               a(k-1)=eps/(dz*dz)
               b(k)=-(ll*ll + kk*kk + 2.*eps/(dz*dz))
               c(k)=eps/(dz*dz)
               rhs(k)=F(i,j,k)
            enddo
               a(P-1)=eps/(dz*dz)
               b(P)=-(ll*ll + kk*kk + eps/(dz*dz))
               rhs(P)=F(i,j,P)
            call tridiag(a,b,c,rhs,P)
C               phibt(i,j)=Fbt(i,j)/(ll*ll + kk*kk)
            do k=1,P
C               phi(i,j,k)=rhs(k)+phibt(i,j)
C               vsp(i,j,k)=im*kk*(rhs(k)+phibt(i,j))
               phi(i,j,k)=rhs(k)
               vsp(i,j,k)=im*kk*rhs(k)
               phi(2*N+2-i,j,k)=conjg(phi(i,j,k))
               vsp(2*N+2-i,j,k)=conjg(vsp(i,j,k))
            enddo            
         enddo !!!! end y-wavenumber loop
      enddo !!!! end x-wavenumber loop

CCCCC do mean flow
      do j=1,M-1
         ll=(j-1.)*pie/Ly
         if (j.eq.1) then
            ubtsp(j)=0.
         else
            ubtsp(j)=-Fybt(j)/(ll*ll)
         endif
            b(1)=-(ll*ll+(eps/(dz*dz)))
            c(1)=eps/(dz*dz)
            rhs(1)=Fy(j,1)
         do k=2,P-1
            a(k-1)=eps/(dz*dz)
            b(k)=-(ll*ll + 2.*eps/(dz*dz))
            c(k)=eps/(dz*dz)
            rhs(k)=Fy(j,k)
         enddo
            a(P-1)=eps/(dz*dz)
            b(P)=-(ll*ll + eps/(dz*dz))
            rhs(P)=Fy(j,P)
         call tridiag(a,b,c,rhs,P)
         do k=1,P
            usp(j,k)=rhs(k) + ubtsp(j)
         enddo
      enddo !!!! end y-wavenumber loop

CCCCC transform to physical space      
      do k=1,P
CCCCC full ifft in x for eddies
         do j=1,M-1
            do i=1,2*N
               fftxr(i)=real(phi(i,j,k))
               fftxi(i)=imag(phi(i,j,k))
               fftxvr(i)=real(vsp(i,j,k))
               fftxvi(i)=imag(vsp(i,j,k))
            enddo
            call sffteu(fftxr,fftxi,2*N,pwr2*2,-1) 
            call sffteu(fftxvr,fftxvi,2*N,pwr2*2,-1) 
            do i=1,2*N,2
               tmp2((i+1)/2,j)=4.*fftxr(i)
               tmp3((i+1)/2,j)=4.*fftxvr(i)
            enddo
         enddo
CCCCC do inverse sin transform in y for mean and eddies
         do i=1,N
            do j=1,M-1
               fftyrp(j)=tmp2(i,j)
               fftyrv(j)=tmp3(i,j)
               fftyrp(j+M-1)=0.
               fftyrv(j+M-1)=0.
            enddo
            call sinft(fftyrp,2*(M-1))                
            call sinft(fftyrv,2*(M-1))                
            do j=1,2*(M-1),2
               psi(i,(j+1)/2,k) = fftyrp(j)/(M-1.)
               v(i,(j+1)/2,k) = fftyrv(j)/(M-1.)
            enddo
               psi(i,M,k) = 0.
               v(i,M,k) = 0.
         enddo
         do j=1,M-1
            fftyr(j)=usp(j,k)
            fftyr(j+M-1)=0.
         enddo         
            call sinft(fftyr,2*(M-1))
         do j=1,2*(M-1),2
            u((j+1)/2,k)=2.*fftyr(j)/(M-1.)
         enddo
            u(M,k)=0.
      enddo

      return
      end

