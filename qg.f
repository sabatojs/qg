      program qg

      implicit none

      integer N,M,P,i,j,k,Tf,tt,pwr2
      parameter(N=128,M=129,P=30)

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

      integer k_in
      real U_in,N_in,H_in,D_in,B_in,f_in

      pie = 4.d0*datan(1.d0)
      im = cmplx(0.,1.)
      dofilter=.true.


      pwr2 = 1
      do while (2**pwr2.lt.N)
         pwr2=pwr2+1
      enddo

      nlin = 1.
      equil = 1.
      eq_eddy = 1.

      open(10,FILE='kx.txt')
      read(10,*) kin
      close(10)
      open(11,FILE='kxchar.txt')
      read(11,*) wavestr
      close(11)

CCCC after read in reset default values to the gcm values      
      H=1.e4
      f=1.e-4
      N2=4.e-4
      B=0.e-11
      d=0.e-18
      U=25.

      k_in = int(kin)
      open(10,FILE='H0.txt')
      open(11,FILE='f0.txt')
      open(12,FILE='N0.txt')
      open(13,FILE='B0.txt')
      open(14,FILE='D0.txt')
      open(15,FILE='U0.txt')
      do i=1,k_in
	 read(10,*) H_in
	 read(11,*) f_in
	 read(12,*) N_in
	 N_in = N_in*N_in
	 read(13,*) B_in
	 read(14,*) D_in
	 read(15,*) U_in
      enddo
      close(10)
      close(11)
      close(12)
      close(13)
      close(14)
      close(15)

C      H=H_in
C      f=f_in
C      N2=N_in
C      B=B_in
C      d=D_in*2.
C      U=U_in

      g=9.81
      TEQ = 305.
      
      Ly = 1.e7
      Lx = 1.e7

 
C      Tf = int(equil)*24*1440 + int(nlin)*4*1440 + 2*1440 + 1
      Tf = int(equil)*24*1440 + int(nlin)*2*1440 + 2*1440 + 1
      wrtflg=0
C      howoften = 36
      howoften = 24
      recno = 1

      eps = f*f/N2
      lam = U/H

      nn = 4.
      nu = nlin*1.e5
      gamma = 0.05

      dx = Lx/N
      dy = Ly/(M-1.)
      dz = H/P

      dt = 900.
C      dt = 3600.

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


CCCCCCCCCCCCCCCCCCCCCC wavenumber loop CCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCC wavenumber loop CCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCC wavenumber loop CCCCCCCCCCCCCCCCCCCCCCCCCCC
      do k0=1,1
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

CCCCCCCCCC eddy initialization CCCCCCCCCCC
      if (nlin.eq.1) then
      if (equil.eq.1.) then
CCCC nonlinear model equilibrium run initialized with 
CCCC random surface temperature perturbation
CCCC details in random_ics amount to a spectrum
CCCC peaked around ktot = (k^2 + l^2 )^1/2 = 2
CCCC with random fourier phases and decaying to zero at
CCCC meridional walls
       write(*,*) 
     $ 'Generating random surface/lid temperature perturbations...'
       call random_ics(q1,N,M,P)
       write(*,*) 'Done!'
      maxq1 = 0.
      do i=1,N
         do j=1,M
            do k=1,P
               if (abs(q1(i,j,k)).gt.maxq1) maxq1=abs(q1(i,j,k))
            enddo
         enddo
      enddo         
      do i=1,N
         do j=1,M
            do k=1,P
               q1(i,j,k)=q1(i,j,k)/maxq1
            enddo
         enddo
      enddo         
      do i=1,N
         do j=1,M
C           read(10,*) q1(i,j,1)
            q1(i,j,1)=10.*(U/Lx)*q1(i,j,1)
            q1(i,j,P)=-q1(i,j,1)
         enddo
      enddo
      else ! non-linear but not equilibrium run
CCCC nonlinear model initialized with linear normal mode
CCCC output to "INIT_FILE" has max(psi)=1 and everything
CCCC else scaled appropriately ==> so scale eddy PV 
CCCC accordingly
C      call ncreadhist(q1,psi1,tmp1,
C     $          tmp2,TIME,1,N,M,P,INIT_FILE,1)
      call ncreadhist(q1,psi1,tmp8,tmp1,
     $          tmp2,y,x,z,
     $          TIME,1,N,M,P,INIT_FILE,1)
      do i=1,N
         do j=1,M
            do k=1,P
               q1(i,j,k)=(5.e-3)*U*Lx*q1(i,j,k)
            enddo
         enddo
      enddo         
      endif ! whether equilibrium or life-cycle run
      else ! linear run
CCCC linear model initialized with cos(k0*x)
      do i=1,N
         do j=1,M
C   	       q1(i,j,1) = - 0.0001*(f*g/N2)*cos(k0*x(i)*2.*pie/Lx)*
CC     $			cos(pie*y(j)/Ly)/TH(1)
C     $		cos(pie*y(j)/Ly)*exp(-64.*y(j)*y(j)/(Ly*Ly))/TH(1)
C            do k=2,P-1
C	       q1(i,j,k) = 0.
C            enddo
C	       q1(i,j,P) = + 0.0001*(f*g/N2)*cos(k0*x(i)*2.*pie/Lx)*
CC     $			cos(pie*y(j)/Ly)/TH(P)
C     $		cos(pie*y(j)/Ly)*exp(-64.*y(j)*y(j)/(Ly*Ly))/TH(P)

               psi0(i,j,1)=cos(2.*pie*x(i)/Lx)*cos(pie*y(j)/Ly)
     $                  *sin(2.*pie*z(1)/H)
               q1(i,j,1) =-pie*pie*(4./(Lx*Lx) + 1./(Ly*Ly))*
     $          cos(2.*pie*x(i)/Lx)*cos(pie*y(j)/Ly)*sin(2.*pie*z(1)/H)
     $          +(sin(2.*pie*z(2)/H)-sin(2.*pie*z(1)/H))
     $          *cos(2.*pie*x(i)/Lx)*cos(pie*y(j)/Ly)*(eps/(dz*dz))
            do k=2,P-1
               psi0(i,j,k)=cos(2.*pie*x(i)/Lx)*cos(pie*y(j)/Ly)
     $                  *sin(2.*pie*z(k)/H)
               q1(i,j,k) = -pie*pie*(4./(Lx*Lx) + 1./(Ly*Ly)
     $          + 4.*f*f/(N2*H*H))*
     $          cos(2.*pie*x(i)/Lx)*cos(pie*y(j)/Ly)*sin(2.*pie*z(k)/H)
            enddo
               psi0(i,j,P)=cos(2.*pie*x(i)/Lx)*cos(pie*y(j)/Ly)
     $                  *sin(2.*pie*z(P)/H)
               q1(i,j,P) =-pie*pie*(4./(Lx*Lx) + 1./(Ly*Ly))*
     $          cos(2.*pie*x(i)/Lx)*cos(pie*y(j)/Ly)*sin(2.*pie*z(P)/H)
     $          -(sin(2.*pie*z(P)/H)-sin(2.*pie*z(P-1)/H))
     $          *cos(2.*pie*x(i)/Lx)*cos(pie*y(j)/Ly)*(eps/(dz*dz))
         enddo      
      enddo
      endif ! linear vs. nonlinear
CCCCC end eddy initialization CCCCC

CCCCC mean flow initialization CCCCC
C damping time in days 6.5*(atan(((-7.5*abs(y) + 10*abs((3/2)*pi./y))))+pi/2)
C (formula above has y in radians from -pi to pi but below y is in meters from 
C -Ly/2 to Ly/2 )
C mean flow is damped back to initial condition on this time scale
C here is where you would change the mean flow ==> U is defined here and
C then Qy must be derived analytically and defined here as well
         mintau=100000.
         ss = Ly/12.
         do j=1,M
            do k=1,P
CCCCCC U = Uz/H * (cos(pi*y/L)^nn)
C               U1(j,k)=lam*z(k)*(cos(pie*y(j)/Ly)**nn)
C               Qy0(j,k)=(nn*pie*pie*lam*z(k)/(Ly*Ly))
C     $		*(cos(pie*y(j)/Ly)**(nn-2.))
C     $		*((cos(pie*y(j)/Ly)**2.)
C     $		-(nn-1.)*(sin(pie*y(j)/Ly)**2.))
CCCCCC U = Uz/H * cos(pi*y/L) * exp( -y^2/(2s^2)), s=L/10
               U1(j,k) = lam*z(k)*cos(pie*y(j)/Ly)*
     $			exp(-y(j)*y(j)/(2.*ss*ss))
               Qy0(j,k)=lam*z(k)*
     $		exp(-y(j)*y(j)/(2.*ss*ss))*
     $		(-(2.*pie*y(j)/(Ly*ss*ss))*sin(pie*y(j)/Ly) 
     $		+ (pie*pie/(Ly*Ly) + 1./(ss*ss) 
     $		- y(j)*y(j)/(ss*ss*ss*ss))*cos(pie*y(j)/Ly))
               tau(j,k) = (atan(((-7.5*abs(y(j)*2.*pie/Ly) 
     $		+ 10.*abs((3./2.)*pie/(y(j)*2.*pie/Ly)))))+pie/2)
               mintau=min(mintau,tau(j,k))
            enddo
         enddo
         maxtau=0.
         do j=1,M
            do k=1,P
               maxtau=max(maxtau,tau(j,k)-mintau)
            enddo
         enddo
C	 kin = sqrt(N2)/1.e-2
         do j=1,M
            do k=1,P
               tau(j,k)=(tau(j,k)-mintau)*(70.-5.)*86400./maxtau  
     $			+ 5.*86400.
            enddo
         enddo

C         tau = tau*sqrt(N2/1.e-4)

CCCC show user what damping times look like CCCC
      if (equil.eq.1.) then
      write(*,*)
      write(*,*) 'DAMPING TIMES:'
      do j=1,M
         write(*,2000) ' ',tau(j,1)/86400.,' at ',y(j)/1000.,' km N'
      enddo
      write(*,*)
      else
      endif

      ! add surface/lid temperature to Qy according to
      ! method of solution following arakawa where approximation
      ! to delta-function PV is integrated in time rather than
      ! explicitly integrating temperature on either boundary
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

      ! given initial q,Qy re-calculate  psi,v,U for 
      ! numerical consistency - also because original
      ! initialization above for linear mode contains
      ! only eddy part but stream calculates total psi,q
      psi1=0.
      v1=0.
      U1=0.
      call stream(psi1,v1,U1,q1,Qy1,
     $		    dy,dz,eps,Lx,Ly,N,M,P,pwr2,nlin,k0)

CCCCC renormalize initial psi,v,q for linear runs
CCCC if we're running linear normal mode calculation
CCCC absolute value of variables becomes very large
CCCC but we may not have converged in growth rates
CCCC before we get Inf or NaN so we need to rescale
         if (nlin.eq.0.) then
         norm1 = 0.
         do i=1,N
            do j=1,M
               do k=1,P
		  norm1 = norm1+psi1(i,j,k)*psi1(i,j,k)
               enddo
            enddo
         enddo
         norm1 = sqrt(norm1)
         do i=1,N
	    do j=1,M
               do k=1,P
	          psi1(i,j,k)=psi1(i,j,k)/norm1
	          q1(i,j,k)=q1(i,j,k)/norm1
	          v1(i,j,k)=v1(i,j,k)/norm1
               enddo
	    enddo
         enddo
         else
         endif ! if linear

      ! for linear runs create netcdf file and output
      ! initial conditions in dimensional units, including
      ! actual temperature not just deviation from standard
      ! atmosphere
      if (nlin.eq.0.) then
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
	       QHIST(i,j,k)=q1(i,j,k)
	       UHIST(j,k)=U1(j,k)
	       QYHIST(j,k)=Qy1(j,k)
            enddo
         enddo
      enddo
      recno = 1
      call ncwritehist_lin(QHIST,PHIST,QYHIST,UHIST,TSHIST,
     $          TTHIST,0.,0.,1,N,M,P,OUTFILELIN,recno)
      else 
      ! nonlinear model output initial condition as well
      call ncinit(OUTFILE,P,M,N,z,y,x,
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
	       QHIST(i,j,k)=q1(i,j,k)
	       UHIST(j,k)=U1(j,k)
	       QYHIST(j,k)=Qy1(j,k)
            enddo
         enddo
      enddo
      recno = 1
      call ncwritehist(QHIST,PHIST,QYHIST,UHIST,TSHIST,
     $          TTHIST,0.,1,N,M,P,OUTFILE,recno)
      endif ! if linear run

CCCC tell user we've started the run
      days = TIME/86400.
      dys = int(days)
      hours = (TIME - dble(dys)*86400.)/3600.
      hrs = int(hours)
      mins = (TIME-dble(dys)*86400.-dble(hrs)*3600.)/60.   
      mint = int(mins)
      secs = (TIME-dble(dys)*86400.-dble(hrs)*3600.-dble(mint)*60.)
      sec = int(secs)
      write(*,3000) ' TIME: ',dys,' days '
     +              ,hrs,' hours '
     +              ,mint,' minutes '
     +              ,sec,' seconds '
      if (nlin.eq.0.) write(*,*) 'GROWTH RATE: N/A for first time step'
      if (nlin.eq.0.) write(*,*) 
     $ 'GROWTH RATE CHANGE (%): N/A for first time step'
      write(*,*)

CCCC start calculating in earnest
CCCC diffuse eddy PV
      do k=1,P
         do i=1,N
            do j=1,M
               tmp1(i,j)=q1(i,j,k)
            enddo
         enddo
         call crank_nicol(tmp1,dx,dy,dt,nu,N,M)
         do i=1,N
            do j=1,M
               q1(i,j,k)=tmp1(i,j)
            enddo
         enddo
CCCC advect eddy PV
            j=1
            tmp7(j)=0.
            tmp9(j)=0.
            tmp8(j,k)=0.
            do i=1,N
               tmp1(i,j) = psi1(i,j,k)
               tmp2(i,j) = q1(i,j,k)
            enddo
         do j=2,M-1
            tmp9(j)=0.
            tmp7(j)=0.
            tmp8(j,k)=0.
            do i=1,N
               tmp1(i,j) = psi1(i,j,k)
               tmp2(i,j) = q1(i,j,k)
               tmp7(j) = ( (i-1.)*tmp7(j) + 
     +			( q1(i,j+1,k)*v1(i,j+1,k)
     +			- q1(i,j-1,k)*v1(i,j-1,k))/(2.*dy)
     +			)/i
               tmp8(j,k) = ( (i-1.)*tmp8(j,k) + 
     +			( q1(i,j+1,k)*v1(i,j+1,k)
     +			+ q1(i,j-1,k)*v1(i,j-1,k)
     +			-2.*q1(i,j,k)*v1(i,j,k) )/(dy*dy)
     +			)/i
               tmp9(j) = ( (i-1.)*tmp9(j) + 
     +			 q1(i,j,k)*v1(i,j,k) 
     +			)/i
            enddo
	 enddo
            j=M
            tmp7(j)=0.
            tmp9(j)=0.
            tmp8(j,k)=0.
            j=1
            tmp7(j)=tmp9(j+1)/dy
            tmp8(j,k)=(tmp7(j+1)-tmp7(j))/dy
            j=M
            tmp7(j)=-tmp9(j-1)/dy
            tmp8(j,k)=(tmp7(j)-tmp7(j-1))/dy
            do i=1,N
               tmp1(i,j) = psi1(i,j,k)
               tmp2(i,j) = q1(i,j,k)
            enddo
         call arakawa(tmp3,tmp1,tmp2,dx,dy,N,M)
         ! if U=0 at meridional walls
         ! then q=q(t=0)=0 for all time
            do i=1,N
               q2(i,1,k)=0.
            enddo
         do j=2,M-1
               i=1
               q2(i,j,k) = q1(i,j,k) + dt*(
     +		nlin*(
     +		-tmp3(i,j)  ! jacobian
     +		+tmp7(j) )  ! pv flux divergence
     +		-(B-d*y(j))*v1(i,j,k)
     +		-Qy1(j,k)*v1(i,j,k)
     +		-U1(j,k)*(q1(2,j,k)-q1(N,j,k))/(2.*dx)
     +          -eq_eddy*q1(i,j,k)/tau(j,k)
     +				 )
            do i=2,N-1
               q2(i,j,k) = q1(i,j,k) + dt*(
     +		nlin*(
     +		-tmp3(i,j)  ! jacobian
     +		+tmp7(j) )  ! pv flux divergence
     +		-(B-d*y(j))*v1(i,j,k)
     +		-Qy1(j,k)*v1(i,j,k)
     +		-U1(j,k)*(q1(i+1,j,k)-q1(i-1,j,k))/(2.*dx)
     +          -eq_eddy*q1(i,j,k)/tau(j,k)
     +				 )
            enddo
               i=N
               q2(i,j,k) = q1(i,j,k) + dt*(
     +		nlin*(
     +		-tmp3(i,j)  ! jacobian
     +		+tmp7(j) )  ! pv flux divergence
     +		-(B-d*y(j))*v1(i,j,k)
     +		-Qy1(j,k)*v1(i,j,k)
     +		-U1(j,k)*(q1(1,j,k)-q1(N-1,j,k))/(2.*dx)
     +          -eq_eddy*q1(i,j,k)/tau(j,k)
     +				 )
	 enddo
         ! if U=0 at meridional walls
         ! then q=q(t=0)=0 for all time
         do i=1,N
            q2(i,M,k)=0.
         enddo
      enddo

CCCC diffuse mean PV gradient
      do k=1,P
         do j=1,M
            tmp7(j)=Qy1(j,k)
         enddo
C         call crank_nicol_mean(tmp7,dy,dt,nu,M)
         do j=1,M
            Qy1(j,k)=tmp7(j)
         enddo

CCCC advect mean PV gradient
            ! if U=0 at meridional walls
            ! then Q=Q(t=0) for all time
	    Qy2(1,k) = Qy1(1,k)
         do j=2,M-1
	    Qy2(j,k) = Qy1(j,k) + dt*(
     +		-nlin*tmp8(j,k)
     +		-equil*(Qy1(j,k)-Qy0(j,k))/tau(j,k)
     +				     )
         enddo
            ! if U=0 at meridional walls
            ! then Q=Q(t=0) for all time
	    Qy2(M,k) = Qy1(M,k)
      enddo

CCCC invert for stream function
      psi2=0.
      v2=0.
      U2=0.
      call stream(psi2,v2,U2,q2,Qy2,
     $		    dy,dz,eps,Lx,Ly,N,M,P,pwr2,nlin,k0)

      ! initialize growth rate to zero
      gr1=0.

CCCCC renormalize initial psi,v,q for linear runs
CCCC if we're running linear normal mode calculation
CCCC absolute value of variables becomes very large
CCCC but we may not have converged in growth rates
CCCC before we get Inf or NaN so we need to rescale
         if (nlin.eq.0.) then
         norm2 = 0.
         do i=1,N
            do j=1,M
               do k=1,P
		  norm2 = norm2+psi2(i,j,k)*psi2(i,j,k)
               enddo
            enddo
         enddo
         norm2 = sqrt(norm2)
         do i=1,N
	    do j=1,M
               do k=1,P
	          psi1(i,j,k)=psi1(i,j,k)/norm2
	          q1(i,j,k)=q1(i,j,k)/norm2
	          v1(i,j,k)=v1(i,j,k)/norm2
	          psi2(i,j,k)=psi2(i,j,k)/norm2
	          q2(i,j,k)=q2(i,j,k)/norm2
	          v2(i,j,k)=v2(i,j,k)/norm2
               enddo
	    enddo
         enddo
         else
         endif ! if linear

CCCCCCCCCCCCCCCCCCCCCC time loop CCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCC time loop CCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCC time loop CCCCCCCCCCCCCCCCCCCCCCCCCCC
      do tt=2,Tf
CCCC diffuse eddy PV
      do k=1,P
         do i=1,N
            do j=1,M
               tmp1(i,j)=q2(i,j,k)
            enddo
         enddo
         call crank_nicol(tmp1,dx,dy,dt,nu,N,M)
         do i=1,N
            do j=1,M
               q2(i,j,k)=tmp1(i,j)
            enddo
         enddo
CCCC advect eddy PV
            j=1
            tmp9(j)=0.
            tmp7(j)=0.
            tmp8(j,k)=0.
            do i=1,N
               tmp1(i,j) = psi2(i,j,k)
               tmp2(i,j) = q2(i,j,k)
            enddo
         do j=2,M-1
            tmp7(j)=0.
            tmp8(j,k)=0.
            do i=1,N
               tmp1(i,j) = psi2(i,j,k)
               tmp2(i,j) = q2(i,j,k)
               tmp7(j) = ( (i-1.)*tmp7(j) + 
     +			( q2(i,j+1,k)*v2(i,j+1,k)
     +			- q2(i,j-1,k)*v2(i,j-1,k))/(2.*dy)
     +			)/i
               tmp8(j,k) = ( (i-1.)*tmp8(j,k) + 
     +			( q2(i,j+1,k)*v2(i,j+1,k)
     +			+ q2(i,j-1,k)*v2(i,j-1,k)
     +			-2.*q2(i,j,k)*v2(i,j,k) )/(dy*dy)
     +			)/i
               tmp9(j) = ( (i-1.)*tmp9(j) + 
     +			 q2(i,j,k)*v2(i,j,k)
     +			)/i
            enddo
	 enddo
            j=M
            tmp7(j)=0.
            tmp9(j)=0.
            tmp8(j,k)=0.
            j=1
            tmp7(j)=tmp9(j+1)/dy
            tmp8(j,k)=(tmp7(j+1)-tmp7(j))/dy
            j=M
            tmp7(j)=-tmp9(j-1)/dy
            tmp8(j,k)=(tmp7(j)-tmp7(j-1))/dy
            do i=1,N
               tmp1(i,j) = psi2(i,j,k)
               tmp2(i,j) = q2(i,j,k)
            enddo
         call arakawa(tmp3,tmp1,tmp2,dx,dy,N,M)
            do i=1,N
	       q3(i,1,k)=0.
            enddo
         do j=2,M-1
               i=1
               q3(i,j,k) = q1(i,j,k) + 2.*dt*(
     +		nlin*(
     +		-tmp3(i,j)  ! jacobian
     +		+tmp7(j) )  ! pv flux divergence
     +		-(B-d*y(j))*v2(i,j,k)
     +		-Qy2(j,k)*v2(i,j,k)
     +		-U2(j,k)*(q2(2,j,k)-q2(N,j,k))/(2.*dx)
     +          -eq_eddy*q2(i,j,k)/tau(j,k)
     +				 )
            do i=2,N-1
               q3(i,j,k) = q1(i,j,k) + 2.*dt*(
     +		nlin*(
     +		-tmp3(i,j)  ! jacobian
     +		+tmp7(j) )  ! pv flux divergence
     +		-(B-d*y(j))*v2(i,j,k)
     +		-Qy2(j,k)*v2(i,j,k)
     +		-U2(j,k)*(q2(i+1,j,k)-q2(i-1,j,k))/(2.*dx)
     +          -eq_eddy*q2(i,j,k)/tau(j,k)
     +				 )
            enddo
               i=N
               q3(i,j,k) = q1(i,j,k) + 2.*dt*(
     +		nlin*(
     +		-tmp3(i,j)  ! jacobian
     +		+tmp7(j) )  ! pv flux divergence
     +		-(B-d*y(j))*v2(i,j,k)
     +		-Qy2(j,k)*v2(i,j,k)
     +		-U2(j,k)*(q2(1,j,k)-q2(N-1,j,k))/(2.*dx)
     +          -eq_eddy*q2(i,j,k)/tau(j,k)
     +				 )
	 enddo
            do i=1,N
	       q3(i,M,k)=0.
            enddo
      enddo

CCCC diffuse mean PV gradient
      do k=1,P
         do j=1,M
            tmp7(j)=Qy2(j,k)
         enddo
         call crank_nicol_mean(tmp7,dy,dt,nu,M)
         do j=1,M
            Qy2(j,k)=tmp7(j)
         enddo

CCCC advect mean PV gradient
	    Qy3(1,k) = Qy1(1,k)
         do j=2,M-1
	    Qy3(j,k) = Qy1(j,k) + 2.*dt*(
     +		-nlin*tmp8(j,k)
     +		-equil*(Qy2(j,k)-Qy0(j,k))/tau(j,k)
     +				     )
         enddo
	    Qy3(M,k) = Qy1(M,k)
      enddo

CCCC invert for stream function
      psi3=0.
      v3=0.
      U3=0.
      call stream(psi3,v3,U3,q3,Qy3,
     $		    dy,dz,eps,Lx,Ly,N,M,P,pwr2,nlin,k0)

CCCCC time filtering to control leap-from computaional mode
      if (dofilter) then
	       q2 = q2 + gamma*(q1+q3-2.*q2)
	       psi2 = psi2 + gamma*(psi1+psi3-2.*psi2)
	       v2 = v2 + gamma*(v1+v3-2.*v2)
	       Qy2 = Qy2 + gamma*(Qy1+Qy3-2.*Qy2)
	       U2 = U2 + gamma*(U1+U3-2.*U2)
      endif

CCCC calculate growth rates
      if (nlin.eq.0.) then
         norm1 = 0.
         norm2 = 0.
         norm3 = 0.
         do i=1,N
            do j=1,M
               do k=1,P
		  norm1 = norm1+psi1(i,j,k)*psi1(i,j,k)
		  norm2 = norm2+psi2(i,j,k)*psi2(i,j,k)
		  norm3 = norm3+psi3(i,j,k)*psi3(i,j,k)
               enddo
            enddo
         enddo
C         norm1 = sqrt(norm1)
C         norm2 = sqrt(norm2)
C         norm3 = sqrt(norm3)
C         gr2 = (norm3-norm1)/(4.*dt*norm2)
         gr2 = (log(norm3)-log(norm1))/(4.*dt)
         tol = (gr2-gr1)*2./(gr1+gr2)
      else
      endif

CCCC update history file
      TIME = dt*(tt-1)
      wrtflg = wrtflg+1

      if (wrtflg.eq.howoften) then
      wrtflg=0
      recno=recno+1
      do k=1,P
            PBAR(1,k) = 0.
	    PBAR(2,k)=PBAR(1,k) - U2(1,k)*dy
         do j=2,M-1
	    PBAR(j+1,k)=PBAR(j,k) - U2(j,k)*dy
         enddo
      enddo
         do i=1,N
            do j=1,M
            TSHIST(i,j)=( ( psi2(i,j,2)-psi2(i,j,1) 
     $		+PBAR(j,2)-PBAR(j,1) ) 
     $		*f/(dz*g)+1.)*TH(1)
            TTHIST(i,j)=( ( psi2(i,j,P)-psi2(i,j,P-1) 
     $		+ PBAR(j,P)-PBAR(j,P-1) ) 
     $		*f/(dz*g)+1.)*TH(P)
               do k=1,P
	          PHIST(i,j,k)=psi2(i,j,k)
	          QHIST(i,j,k)=q2(i,j,k)
	          UHIST(j,k)=U2(j,k)
	          QYHIST(j,k)=Qy2(j,k)
               enddo
            enddo
         enddo

      if (nlin.eq.0.) then
      call ncwritehist_lin(QHIST,PHIST,QYHIST,UHIST,TSHIST,
     $          TTHIST,gr2,TIME,tt,N,M,P,OUTFILELIN,recno)
      else
      call ncwritehist(QHIST,PHIST,QYHIST,UHIST,TSHIST,
     $          TTHIST,TIME,tt,N,M,P,OUTFILE,recno)
      endif ! if linear run

      days = TIME/86400.
      dys = int(days)
      hours = (TIME - dble(dys)*86400.)/3600.
      hrs = int(hours)
      mins = (TIME-dble(dys)*86400.-dble(hrs)*3600.)/60.   
      mint = int(mins)
      secs = (TIME-dble(dys)*86400.-dble(hrs)*3600.-dble(mint)*60.)
      sec = int(secs)
         
      write(*,3000) ' TIME: ',dys,' days '
     +              ,hrs,' hours '
     +              ,mint,' minutes '
     +              ,sec,' seconds '
      if (nlin.eq.0.) write(*,*) 'GROWTH RATE: ',gr2
      if (nlin.eq.0.) write(*,*) 'GROWTH RATE CHANGE (%): ',tol*100.
      write(*,*)

      endif ! if it's time to output

CCCC update variables
	       q1 = q2
	       q2 = q3
	       q3 = 0.
	       psi1 = psi2
	       psi2 = psi3
	       psi3 = 0.
	       v1 = v2
	       v2 = v3
	       v3 = 0.
	       Qy1 = Qy2
	       Qy2 = Qy3
	       Qy3 = 0.
	       U1 = U2
	       U2 = U3
	       U3 = 0.
               gr1=gr2

CCCC if we're running linear normal mode calculation
CCCC absolute value of variables becomes very large
CCCC but we may not have converged in growth rates
CCCC before we get Inf or NaN so we need to rescale
         if (nlin.eq.0.) then
C         maxpsi=0.
C         do i=1,N
C	    do j=1,M
C               do k=1,P
C	          maxpsi=max(abs(psi1(i,j,k)),maxpsi)
C               enddo
C	    enddo
C         enddo
          maxpsi = norm1
C         if (maxpsi.ge.1.e10) then
         do i=1,N
	    do j=1,M
               do k=1,P
	          psi1(i,j,k)=psi1(i,j,k)/maxpsi
	          q1(i,j,k)=q1(i,j,k)/maxpsi
	          v1(i,j,k)=v1(i,j,k)/maxpsi
	          psi2(i,j,k)=psi2(i,j,k)/maxpsi
	          q2(i,j,k)=q2(i,j,k)/maxpsi
	          v2(i,j,k)=v2(i,j,k)/maxpsi
               enddo
	    enddo
         enddo
C         else
C         endif ! if psi got too big
         else
         endif ! if linear

      enddo ! time loop
CCCCCCCCCCCCCCCCCCCCCC time loop CCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCC time loop CCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCC time loop CCCCCCCCCCCCCCCCCCCCCCCCCCC


CCCCCC rescale and write out variables at end
CCCCCC of linear run to initialization file to
CCCCCC be read by nonlinear calculation
      if (nlin.eq.0) then
      call ncinit(INIT_FILE,P,M,N,z,y,x,
     $                  TH,N2,H,B,d,f,g,U,Lx,Ly,gr2,k0*1.)
         maxpsi=0.
         do i=1,N
	    do j=1,M
               do k=1,P
	          maxpsi=max(psi2(i,j,k),maxpsi)
               enddo
	    enddo
         enddo
         do i=1,N
	    do j=1,M
               do k=1,P
	          psi2(i,j,k)=psi2(i,j,k)/maxpsi
	          q2(i,j,k)=q2(i,j,k)/maxpsi
               enddo
	    enddo
         enddo
      do i=1,N
         do j=1,M
            TSHIST(i,j)=( ( psi2(i,j,2)-psi2(i,j,1) )
     $		*f/(dz*g)+1.)*TH(1)
            TTHIST(i,j)=( ( psi2(i,j,P)-psi2(i,j,P-1) )
     $		*f/(dz*g)+1.)*TH(P)
            do k=1,P
	       PHIST(i,j,k)=psi2(i,j,k)
	       QHIST(i,j,k)=q2(i,j,k)
	       UHIST(j,k)=U2(j,k)
	       QYHIST(j,k)=Qy2(j,k)
            enddo
         enddo
      enddo
      recno = 1
      call ncwritehist(QHIST,PHIST,QYHIST,UHIST,TSHIST,
     $          TTHIST,0.,1,N,M,P,INIT_FILE,recno)
      else
      endif

CCCCCC reset all variables to zero for next zonal wavenumber
      do j=1,M
         do k=1,P
            do i=1,N
	       psi1(i,j,k)=0.
	       q1(i,j,k)=0.
	       v1(i,j,k)=0.
	       psi2(i,j,k)=0.
	       q2(i,j,k)=0.
	       v2(i,j,k)=0.
	       psi3(i,j,k)=0.
	       q3(i,j,k)=0.
	       v3(i,j,k)=0.
            enddo
	       U1(j,k)=0.
	       Qy1(j,k)=0.
	       U2(j,k)=0.
	       Qy2(j,k)=0.
	       U3(j,k)=0.
	       Qy3(j,k)=0.
	 enddo
      enddo
         gr1=0.
         gr2=0.

      enddo ! wavenumber loop
CCCCCCCCCCCCCCCCCCCCCC wavenumber loop CCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCC wavenumber loop CCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCC wavenumber loop CCCCCCCCCCCCCCCCCCCCCCCCCCC

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

