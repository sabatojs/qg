      subroutine nrg_prod(N,M,T,eps)

      implicit none

CCCC input CCCC
      integer N,M,T

CCCC output CCCC
      real eps(N,M,T)
      
CCCC local variables CCCC
      integer kk,ll,ww
      integer Ktot,Wtot
      real k(N),l(M),w(T)
      real norm,W0,K0
      integer pwr2

      integer randseed
      real randnum
      real maxeps
      real ran1
 
      real rfftw(T),ifftw(T)
      real rfftx(N),ifftx(N)
      real rffty(M),iffty(M)

      complex im
      real pie

      complex spec(N,M,T),specxy(N,M,T),specx(N,M,T)
      complex spec1(N,M,T)
      real eps1(N,M,T)

      write(*,*) 'Enter a random number (integer):'
      read(*,*) randseed
C      randseed = -23

      im = cmplx(0.,1.)
      pie = 4.*atan(1.)

      k(1) = 0.
      do kk = 2,N/2+1
         k(kk) = (kk-1.)
      enddo
      do kk = N/2+2,N
         k(kk) = -N + kk - 1.
      enddo

      l(1) = 0.
      do ll = 2,M/2+1
         l(ll) = (ll-1.)
      enddo
      do ll = M/2+2,M
         l(ll) = -M + ll - 1.
      enddo

      w(1) = 0.
      do ww = 2,T/2+1
         w(ww) = (ww-1.)
      enddo
      do ww = T/2+2,T
         w(ww) = -T + ww - 1.
      enddo


C      write(*,*) 'SEED FOR RANDOM NUMBER GENERATOR:  ',
C     $			randseed
C      randnum = ran1(randseed)
C      write(*,*) 'UNIFORM DEVIATE:  ',randnum
C      write(*,*) 'SEED FOR RANDOM NUMBER GENERATOR:  ',
C     $			randseed
      
      Ktot = N/4
      Wtot = Ktot*T/N !!!! assumes N=M

      do kk=1,N
         do ll=1,M
            do ww=1,T
	       randnum=ran1(randseed)
               K0 = sqrt(k(kk)*k(kk)+l(ll)*l(ll))
               W0 = sqrt(w(ww)*w(ww))
	       spec1(kk,ll,ww) =
     $		exp( cmplx(-0.001*(K0-Ktot)*(K0-Ktot)
     $			    /(2.*pie/(N*N)+2.*pie/(M*M))
     $			   -0.01*(W0-Wtot)*(W0-Wtot)
     $			   ,0.) )
               spec(kk,ll,ww) = spec1(kk,ll,ww)*cexp(
     $			cmplx(0.,2.*pie*randnum))
            enddo
	 enddo
      enddo

      pwr2 = 1
      do while (2**pwr2.lt.T)
         pwr2=pwr2+1
      enddo
      do kk=1,N
         do ll=1,M
            do ww=1,T
               rfftw(ww) = real(spec(kk,ll,ww))
               ifftw(ww) = imag(spec(kk,ll,ww))
            enddo
            call sffteu( rfftw, ifftw, T, pwr2, -1 )
            do ww=1,T
               specxy(kk,ll,ww)=cmplx(rfftw(ww),ifftw(ww))
            enddo
         enddo
      enddo

      pwr2 = 1
      do while (2**pwr2.lt.M)
         pwr2=pwr2+1
      enddo
      do kk=1,N
         do ww=1,T
            do ll=1,M
               rffty(ll) = real(specxy(kk,ll,ww))
               iffty(ll) = imag(specxy(kk,ll,ww))
            enddo
            call sffteu( rffty, iffty, M, pwr2, -1 )
	    do ll=1,M   
               specx(kk,ll,ww)=cmplx(rffty(ll),iffty(ll))
      	    enddo
         enddo
      enddo

      pwr2 = 1
      do while (2**pwr2.lt.N)
         pwr2=pwr2+1
      enddo
      do ll=1,M
         do ww=1,T
            do kk=1,N
               rfftx(kk) = real(specx(kk,ll,ww))
               ifftx(kk) = imag(specx(kk,ll,ww))
            enddo
            call sffteu( rfftx, ifftx, N, pwr2, -1 )
	    do kk=1,N   
               eps1(kk,ll,ww)=rfftx(kk)*(N*M*T)/2.
      	    enddo
         enddo
      enddo

      maxeps = 0.
      do kk=1,N
         do ll=1,M
            do ww=1,T
               norm = abs(eps1(kk,ll,ww))
	       if (norm.gt.maxeps) then
                  maxeps=norm
               else
               endif
            enddo
         enddo
      enddo

      write(*,*) 'MAXIMUM NORM FOR EPSILON:  ',maxeps

      do kk=1,N
         do ll=1,M
            do ww=1,T
               eps(kk,ll,ww)=eps1(kk,ll,ww)/(2.*pie*maxeps)
            enddo
         enddo
      enddo



      return
      end
