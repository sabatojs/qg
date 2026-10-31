      subroutine elliptic(psi,v,u,psi0,u0,q,qy,u1,u2,
     $			dy,dz,eps,L,N,M,P,pwr2,nlin,k0)

CCCC input arguments
      integer N,M,P,pwr2,k0
      real psi0(N,M,P),q(N,M,P),u0(M,P),qy(M,P)
      real u1(P),u2(P)
      real dy,dz,eps,L,nlin

CCCC output arguments
      real psi(N,M,P),v(N,M,P),u(M,P)

CCCC local workspace
      integer i,j,k,flag,maxiters
      real tol,tol1,kk,hh,h2
      complex F(N,M,P),phi(2*N,M,P),vsp(2*N,M,P)
      complex phi0(N,M,P),phi1(N,M,P)
      complex tmp1(N),tmp2(N)
      real fftr1(N),ffti1(N),fftr2(N),ffti2(N)
      real fftr3(N),ffti3(N),fftr4(N),ffti4(N)
      real fftr(2*N),ffti(2*N),fftvr(2*N),fftvi(2*N)
      real NORM,NORM1,RNORM,RNORM1
      complex im
      real pie
      integer k1,k2

      pie = 4.*atan(1.)
      im = cmplx(0.,1.)


      if (nlin.eq.1.) then
         k1=1
         k2=N/2
      else
         k1=k0
         k2=k0+1
      endif

CCCC get fourier transforms
      do i=1,M
         do j=1,P
            do k=1,N
               fftr1(k)=q(k,i,j)
               ffti1(k)=0.0
               fftr2(k)=psi0(k,i,j)
               ffti2(k)=0.0
               phi(2*k-1,i,j)=0.
               phi(2*k,i,j)=0.
               vsp(2*k-1,i,j)=0.
               vsp(2*k,i,j)=0.
            enddo
               call sffteu(fftr1,ffti1,N,pwr2,1)       
               call sffteu(fftr2,ffti2,N,pwr2,1)       
            do k=2,N/2
               F(k,i,j)=cmplx(fftr1(k),ffti1(k))
               phi0(k,i,j)=cmplx(fftr2(k),ffti2(k))
            enddo
               F(1,i,j)=-qy(i,j)
               phi0(1,i,j)=u0(i,j)
         enddo
      enddo

      maxiters = 0

CCCCC do eddies
      do k=1,k2,k1
         maxiters=0
         kk=(k-1.)*2.*pie/L
         hh = 2.*dy*dy*eps/(dz*dz) + 2. + dy*dy*kk*kk
         h2 = dy*dy*eps/(dz*dz) + 2. + dy*dy*kk*kk
         flag=1
      do while (flag.eq.1) !!!! jacobi iteration
         !!!! interior
         do i=2,M-1
            !!!! surface
            j=1
               phi1(k,i,j) = ( -dy*dy*F(k,i,j)         +          (
     $	       phi0(k,i-1,j)+phi0(k,i+1,j)
     $	      +phi0(k,i,j+1)*(dy*dy*eps/(dz*dz))                  )
     $			     )/h2
            do j=2,P-1
               phi1(k,i,j) = ( -dy*dy*F(k,i,j)         +          (
     $	       phi0(k,i-1,j)+phi0(k,i+1,j)
     $	      +(phi0(k,i,j-1)+phi0(k,i,j+1))*(dy*dy*eps/(dz*dz))  )
     $			     )/hh
            enddo
            !!!! lid
            j=P
               phi1(k,i,j) = ( -dy*dy*F(k,i,j)         +          (
     $	       phi0(k,i-1,j)+phi0(k,i+1,j)
     $	      +phi0(k,i,j-1)*(dy*dy*eps/(dz*dz))                  )
     $			     )/h2
         enddo
         if (k.eq.1) then
         !!! equatorward wall
         i=1
            do j=1,P
               phi1(k,1,j) = u1(j)
            enddo
         !!! polewardward wall
         i=M
            do j=1,P
               phi1(k,M,j) = u2(j)
            enddo
         else
         !!! equatorward wall
         i=1
            do j=1,P
               phi1(k,1,j) = 0.
            enddo
         !!! polewardward wall
         i=M
            do j=1,P
               phi1(k,M,j) = 0.
            enddo
         endif

            RNORM=0.
            NORM=0.
      do i=2,M-1
         j=1
            NORM1 = phi1(k,i,j)*conjg(phi1(k,i,j))
            RNORM1 = 
     $		      (-F(k,i,j)-kk*kk*phi1(k,i,j)+(
     $	    (phi1(k,i+1,j)+phi1(k,i-1,j)-2.*phi1(k,i,j))/(dy*dy)
     $	   +(phi1(k,i,j+1)-phi1(k,i,j))*eps/(dz*dz)
     $		      )                            )
     $			*conjg
     $		      (-F(k,i,j)-kk*kk*phi1(k,i,j)+(
     $	    (phi1(k,i+1,j)+phi1(k,i-1,j)-2.*phi1(k,i,j))/(dy*dy)
     $	   +(phi1(k,i,j+1)-phi1(k,i,j))*eps/(dz*dz)
     $		      )                            )
	    RNORM= max(RNORM, RNORM1)
            NORM = max(NORM, NORM1)
         do j=2,P-1
            NORM1 = phi1(k,i,j)*conjg(phi1(k,i,j))
            RNORM1 = 
     $		      (-F(k,i,j)-kk*kk*phi1(k,i,j)+(
     $	    (phi1(k,i+1,j)+phi1(k,i-1,j)-2.*phi1(k,i,j))/(dy*dy)
     $	   +(phi1(k,i,j+1)+phi1(k,i,j-1)-2.*phi1(k,i,j))*eps/(dz*dz)
     $		      )                            )
     $			*conjg
     $		      (-F(k,i,j)-kk*kk*phi1(k,i,j)+(
     $	    (phi1(k,i+1,j)+phi1(k,i-1,j)-2.*phi1(k,i,j))/(dy*dy)
     $	   +(phi1(k,i,j+1)+phi1(k,i,j-1)-2.*phi1(k,i,j))*eps/(dz*dz)
     $		      )                            )
	    RNORM= max(RNORM, RNORM1)
            NORM = max(NORM, NORM1)
         enddo
         j=P
            NORM1 = phi1(k,i,j)*conjg(phi1(k,i,j))
            RNORM1 = 
     $		      (-F(k,i,j)-kk*kk*phi1(k,i,j)+(
     $	    (phi1(k,i+1,j)+phi1(k,i-1,j)-2.*phi1(k,i,j))/(dy*dy)
     $	   -(phi1(k,i,j)-phi1(k,i,j-1))*eps/(dz*dz)
     $		      )                            )
     $			*conjg
     $		      (-F(k,i,j)-kk*kk*phi1(k,i,j)+(
     $	    (phi1(k,i+1,j)+phi1(k,i-1,j)-2.*phi1(k,i,j))/(dy*dy)
     $	   -(phi1(k,i,j)-phi1(k,i,j-1))*eps/(dz*dz)
     $		      )                            )
	    RNORM= max(RNORM, RNORM1)
            NORM = max(NORM, NORM1)
      enddo

      tol=RNORM/NORM
      if (k.eq.1) tol1=1.e-30
      if (k.ne.1) tol1=1.e-20
      if (tol.ge.tol1) then
      do i=1,M
         do j=1,P
	    phi0(k,i,j)=phi1(k,i,j)
         enddo
      enddo
      maxiters=maxiters+1
      else
      !!! make sure psi is real
      do i=1,M
         do j=1,P
            phi(k,i,j)=phi1(k,i,j)
            if (k.gt.1) phi(2*N+2-k,i,j)=conjg(phi1(k,i,j))
            vsp(k,i,j)=im*kk*phi1(k,i,j)
            if (k.gt.1) vsp(2*N+2-k,i,j)=conjg(vsp(k,i,j))
         enddo
      enddo
      flag=0
      endif
      enddo !!!! jacobi iteration

C      write(*,*) 'k = ',k-1,'   iterations: ',maxiters
      enddo !!!! end wavenumber loop

CCCCC transform to physical space      
      do i=1,M
         do j=1,P
            u(i,j)=phi1(1,i,j)
            phi(1,i,j)=0.
            do k=1,2*N
               fftr(k)=real(phi(k,i,j))
               ffti(k)=imag(phi(k,i,j))
               fftvr(k)=real(vsp(k,i,j))
               fftvi(k)=imag(vsp(k,i,j))
            enddo
               call sffteu(fftr,ffti,2*N,2*pwr2,-1)       
               call sffteu(fftvr,fftvi,2*N,2*pwr2,-1)       
            do k=1,2*N,2
               psi((k+1)/2,i,j) = 2.*fftr(k)
               v((k+1)/2,i,j) = 2.*fftvr(k)
            enddo
         enddo
      enddo

      return
      end

