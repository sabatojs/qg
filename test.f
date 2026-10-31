      program qg

      implicit none

      integer N,M,P,i,j,k,Tf,tt,pwr2
      parameter(N=32,M=33,P=20)

      real dx,dy,dz,dt
      real x(N),y(M),z(P)
      real psi1(N,M,P),q1(N,M,P),v1(N,M,P)
      real psi2(N,M,P),q2(N,M,P),v2(N,M,P)
      real psi3(N,M,P),q3(N,M,P),v3(N,M,P)
      real psi0(N,M,P)
      real U1(M,P),U2(M,P),U3(M,P),U0
      real Qy1(M,P),Qy2(M,P),Qy3(M,P)

      real QHIST(N,M,P),PHIST(N,M,P)
      real QYHIST(M,P),UHIST(M,P)
      real TSHIST(N,M),TTHIST(N,M)
      real TH(P),PBAR(M,P)

      real B,N2,f,g,H,Ld,Lx,Ly,U,TEQ,nu,d
      real TIME
      real gamma
      real gr1,gr2,tol,norm1,norm2,norm3

      real tau(M,P),Qy0(M,P)

      real maxpsi,maxtau,mintau,nn

      real mins,hours,days,secs
      integer mint,hrs,dys,sec

      integer k0
      real nlin,equil,eq_eddy
      real lam,eps
 
      integer wrtflg,howoften,recno

      real tmp1(N,M)
      real tmp2(N,M)
      real tmp3(N,M)
      real tmp4(N,M)
      real tmp5(N,M)
      real tmp6(N,M)
      real tmp7(M)
      real tmp8(M,P)
      real tmp9(M)

      real tmp11(M-1)
      real tmp12(M-1)

      real maxq1,ss

      real pie
      complex im
      character(11) OUTFILE
      character(15) OUTFILELIN
      character(13) INIT_FILE
      character(4) wavestr

      logical dofilter

      real kin

      pie = 4.d0*datan(1.d0)
      im = cmplx(0.,1.)
      dofilter=.true.


      pwr2 = 1
      do while (2**pwr2.lt.N)
         pwr2=pwr2+1
      enddo

      nlin = 0.
      equil = 0.
      eq_eddy = 0.

      H=1.e4
      f=1.e-4
      N2=1.e-4
      g=9.81
      TEQ = 305.
      B=0.e-11
      d=0.e-18
      Ly = 1.e7
      open(10,FILE='kx.txt')
      read(10,*) kin
      close(10)
      open(11,FILE='kxchar.txt')
      read(11,*) wavestr
      close(11)
C      Lx = 5.e7/kin
C      kin = 1.
C      if (kin.eq.1.) Lx = (2.e7)
C      if (kin.eq.1.) wavestr = 'equi'
C      if (kin.eq.2.) Lx = (1.5e7)/2.
C      if (kin.eq.2.) wavestr = '7500'
C      if (kin.eq.3.) Lx = (2.e7)/3.
C      if (kin.eq.3.) wavestr = '6667'
C      if (kin.eq.4.) Lx = (1.e7)/2.
C      if (kin.eq.4.) wavestr = '5000'
C      if (kin.eq.5.) Lx = (1.e7)/3.
C      if (kin.eq.5.) wavestr = '3333'

      Lx = (1.e6)*(21.-kin)/2.

 
      Tf = int(equil)*24*1440 + int(nlin)*4*1440 + 1*1440 + 1
C      Tf = int(equil)*24*1440 + int(nlin)*4*1440 + 4*1440 + 1
      wrtflg=0
C      howoften = 24
      howoften = 6
      recno = 1

      eps = f*f/N2
      lam = 2.e-3
      U=H*lam

      nn = 8.
      nu = nlin*1.e5
      gamma = 0.05

      dx = Lx/N
      dy = Ly/(M-1.)
      dz = H/P

      dt = 3600.

      write(*,*)
      write(*,*)
      write(*,*) 'MODEL PARAMETERS'
      write(*,1110) 'WIND SPEED:          ',U,'       m/s'
      write(*,1200) 'BETA:                ',B,'  1/ms'
      write(*,1200) 'DELTA:               ',d,'  1/m2s'
      write(*,1310) 'CHANNEL LENGTH:      ',dx*N/1000.,'    km'
      write(*,1300) 'CHANNEL WIDTH:       ',dy*(M-1)/1000.,'     km'
      write(*,1300) 'ROSSBY RADIUS:       ',sqrt(N2)*H/(1000.*f)
     $ ,'     km'
      write(*,1400) 'BV FREQUENCY:        ',sqrt(N2),'  1/s'
      write(*,1500) 'DIFFUSIVITY:         ',nu,'  m2/s'
      write(*,1120) 'TIME STEP:           ',dt,'     s'
      write(*,1120) 'MAX TIME STEP (CFL): ',
     $dx/max(U,B*Lx*Lx/(4.*pie*pie)),'    s'
      write(*,1130) 'GRID SPACING: '
      write(*,2100) 
     $	'            DX x DY  ',dx/1000.,' km  x ',dy/1000,' km'
      write(*,1140) '                 DZ  ',dz/1000.,' km'
      write(*,*)
      write(*,*)

      write(*,*)
      write(*,*) 'STANDARD ATMOSPHERE TEMPERATURES:'

CCCCC grid CCCCC
      do i=1,N
         x(i) = dx*(i-1.) - Lx/2.
      enddo
      do j=1,M
         y(j) = dy*(j-1.) - Ly/2.
      enddo
      do k=1,P
         z(k)=(k- 1./2.)*dz
         TH(k) = TEQ*exp(N2*z(k)/g)
         write(*,1000) ' ',TH(k),' K at ',z(k)/1000.,' km'
      enddo


         psi1=0.
         q1=0.
         v1=0.
         Qy1=0.
         U1=0.
         psi2=0.
         q2=0.
         v2=0.
         Qy2=0.
         U2=0.
         psi3=0.
         q3=0.
         v3=0.
         Qy3=0.
         U3=0.

         Qy0=0.
         tau=0.

         tmp1=0.
         tmp2=0.
         tmp3=0.
         tmp4=0.
         tmp5=0.
         tmp6=0.
         tmp7=0.
         tmp8=0.
         tmp9=0.
         tmp11=0.
         tmp12=0.

         INIT_FILE='init_' // wavestr // '.cdf'
         OUTFILE='qg_' // wavestr // '.cdf'
         OUTFILELIN='qg_lin_' // wavestr // '.cdf'
         write(*,*) 
         write(*,*) 'ZONAL WAVENUMBER: ',k0
         write(*,*) 
         if (nlin.eq.0.) write(*,*) 'THIS IS A LINEAR RUN'
         if (nlin.eq.1.) write(*,*) 'THIS IS A NONLINEAR RUN'
         write(*,*) 

      k0 = 1
CCCC linear model initialized with cos(k0*x)
      if (0.eq.1) then
      do i=1,N
         do j=1,M
   	       q1(i,j,1) = - 0.0001*(f*g/N2)*cos(k0*x(i)*2.*pie/Lx)*
     $			cos(pie*y(j)/Ly)/TH(1)
C     $		cos(pie*y(j)/Ly)*exp(-64.*y(j)*y(j)/(Ly*Ly))/TH(1)
            do k=2,P-1
	       q1(i,j,k) = 0.
            enddo
	       q1(i,j,P) = + 0.0001*(f*g/N2)*cos(k0*x(i)*2.*pie/Lx)*
     $			cos(pie*y(j)/Ly)/TH(1)
C     $		cos(pie*y(j)/Ly)*exp(-64.*y(j)*y(j)/(Ly*Ly))/TH(P)
         enddo      
      enddo
      else
      do i=1,N
         do j=1,M
               psi0(i,j,1)=cos(2.*pie*x(i)/Lx)*cos(pie*y(j)/Ly)
     $			*sin(2.*pie*z(1)/H)
   	       q1(i,j,1) =-pie*pie*(4./(Lx*Lx) + 1./(Ly*Ly))*
     $		cos(2.*pie*x(i)/Lx)*cos(pie*y(j)/Ly)*sin(2.*pie*z(1)/H)
     $		+(sin(2.*pie*z(2)/H)-sin(2.*pie*z(1)/H))
     $		*cos(2.*pie*x(i)/Lx)*cos(pie*y(j)/Ly)*(eps/(dz*dz))
            do k=2,P-1
               psi0(i,j,k)=cos(2.*pie*x(i)/Lx)*cos(pie*y(j)/Ly)
     $			*sin(2.*pie*z(k)/H)
	       q1(i,j,k) = -pie*pie*(4./(Lx*Lx) + 1./(Ly*Ly) 
     $		+ 4.*f*f/(N2*H*H))*
     $		cos(2.*pie*x(i)/Lx)*cos(pie*y(j)/Ly)*sin(2.*pie*z(k)/H)
            enddo
               psi0(i,j,P)=cos(2.*pie*x(i)/Lx)*cos(pie*y(j)/Ly)
     $			*sin(2.*pie*z(P)/H)
   	       q1(i,j,P) =-pie*pie*(4./(Lx*Lx) + 1./(Ly*Ly))*
     $		cos(2.*pie*x(i)/Lx)*cos(pie*y(j)/Ly)*sin(2.*pie*z(P)/H)
     $		-(sin(2.*pie*z(P)/H)-sin(2.*pie*z(P-1)/H))
     $		*cos(2.*pie*x(i)/Lx)*cos(pie*y(j)/Ly)*(eps/(dz*dz))
         enddo      
      enddo
      endif
CCCCC end eddy initialization CCCCC

CCCCC mean flow initialization CCCCC
C damping time in days 6.5*(atan(((-7.5*abs(y) + 10*abs((3/2)*pi./y))))+pi/2)
C (formula above has y in radians from -pi to pi but below y is in meters from 
C -Ly/2 to Ly/2 )
C mean flow is damped back to initial condition on this time scale
C here is where you would change the mean flow ==> U is defined here and
C then Qy must be derived analytically and defined here as well
         do j=1,M
            do k=1,P
CCCCCC U = Uz/H * (cos(pi*y/L)^nn)
               U1(j,k)=lam*z(k)*(cos(pie*y(j)/Ly)**nn)
               Qy0(j,k)=(nn*pie*pie*lam*z(k)/(Ly*Ly))
     $		*(cos(pie*y(j)/Ly)**(nn-2.))
     $		*((cos(pie*y(j)/Ly)**2.)
     $		-(nn-1.)*(sin(pie*y(j)/Ly)**2.))
            enddo
         enddo


      ! add surface/lid temperature to Qy according to
      ! method of solution following arakawa where approximation
      ! to delta-function PV is integrated in time rather than
      ! explicitly integrating temperature on the boundaries
      do j=1,M
	    Qy1(j,1) = Qy0(j,1)
     $		- (f*f/N2)*(U1(j,2)-U1(j,1))/(dz*dz)
            Qy0(j,1)=Qy1(j,1)
         do k=2,P-1
	    Qy1(j,k) = Qy0(j,k)
            Qy0(j,k) = Qy1(j,k)
         enddo
	    Qy1(j,P) = Qy0(j,P)
     $		+ (f*f/N2)*(U1(j,P)-U1(j,P-1))/(dz*dz)
            Qy0(j,P)=Qy1(j,P)
      enddo

      psi1=0.
      v1=0.
      U1=0.
      call stream(psi1,v1,U1,q1,Qy1,
     $		    dy,dz,eps,Lx,Ly,N,M,P,pwr2,nlin,k0)

      call ncinit_lin(OUTFILELIN,P,M,N,z,y,x,
     $                  TH,N2,H,B,d,f,g,U,Lx,Ly,0.,k0*1.)
      do k=1,P
            PBAR(1,k) = 0.
	    PBAR(2,k)=PBAR(1,k) - U1(1,k)*dy
         do j=2,M-1
	    PBAR(j+1,k)=PBAR(j,k) - U1(j,k)*dy
         enddo
      enddo
      do i=1,N
         do j=1,M
            TSHIST(i,j)=( ( psi1(i,j,2)-psi1(i,j,1)  
     $		+PBAR(j,2)-PBAR(j,1) )
     $		*f/(dz*g)+1.)*TH(1)
            TTHIST(i,j)=( ( psi1(i,j,P)-psi1(i,j,P-1) 
     $		+PBAR(j,P)-PBAR(j,P-1) )
     $		*f/(dz*g)+1.)*TH(P)
            do k=1,P
	       PHIST(i,j,k)=psi1(i,j,k)
	       QHIST(i,j,k)=psi0(i,j,k)
	       UHIST(j,k)=U1(j,k)
	       QYHIST(j,k)=Qy1(j,k)
            enddo
         enddo
      enddo
      recno = 1
      call ncwritehist_lin(QHIST,PHIST,QYHIST,UHIST,TSHIST,
     $          TTHIST,0.,0.,1,N,M,P,OUTFILELIN,recno)



 1000 format (A,F5.1,A,F5.2,A)
 1110 format (A,F5.1,A)
 1120 format (A,F8.1,A)
 1130 format (A)
 1140 format (A,F7.3,A)
 1100 format (A,F15.10,A)
 1200 format (A,ES10.2,A)
 1300 format (A,F7.1,A)
 1310 format (A,F8.1,A)
 1400 format (A,ES10.2,A)
 1500 format (A,ES10.2,A)
 2000 format (A,F5.2,A,F8.1,A)
 2100 format (A,F6.1,A,F5.1,A)
 3000 format (A,I4,A,I2,A,I2,A,I2,A)
 4000 format (A,F5.0,A)
 9999 format (A,F14.4,A,F12.4)

      return
      end

