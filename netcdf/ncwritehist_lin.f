C     This is part of the netCDF package.
C     Copyright 2006 University Corporation for Atmospheric Research/Unidata.
C     See COPYRIGHT file for conditions of use.

      subroutine ncwritehist_lin(Q,PSI,QY,U,TS,
     $		TT,gr,tim,timestep,N,M,P,FILE_NAME,recno)
      implicit none
      include 'netcdf.inc'

C     This is the name of the data file we will create.
      character*(*) FILE_NAME
      integer ncid

      integer NDIMS
      parameter (NDIMS = 4)
      integer M, N, P

C     The start and count arrays will tell the netCDF library where to
C     write our data.
      integer start(NDIMS), count(NDIMS)
      integer start2(NDIMS-1), count2(NDIMS-1)
      integer start3(NDIMS-1), count3(NDIMS-1)

C      integer t_varid
      integer q_varid
      integer psi_varid
      integer u_varid
      integer qy_varid
      integer tt_varid
      integer ts_varid
      integer tim_varid
      integer gr_varid

C     Program variables to hold the data we will write out. We will only
C     need enough space to hold one timestep of data; one record.
C      real T(N, M, P)
      real Q(N, M, P)
      real U(M, P)
      real QY(M, P)
      real PSI(N, M, P)
      real TS(N,M)
      real TT(N,M)
      real tim,gr

      real*8 mins,hours,days,secs
      integer min,hrs,dys,sec


C     
      integer recno,timestep

C     Error handling.
      integer retval

C     open the file. 
      retval = nf_open(FILE_NAME, nf_write, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     These settings tell netcdf to write one timestep of data. (The
C     setting of start(3) from the time step input tells netCDF which
C     timestep to write.)
      count(1) = N
      count(2) = M
      count(3) = P
      count(4) = 1
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = recno

      count2(1) = N
      count2(2) = M
      count2(3) = 1
      start2(1) = 1
      start2(2) = 1
      start2(3) = recno

      count3(1) = M
      count3(2) = P
      count3(3) = 1
      start3(1) = 1
      start3(2) = 1
      start3(3) = recno

CCCCC lat-lon-lev variables
C         retval = nf_inq_varid(ncid, 'T', t_varid)
C         if (retval .ne. nf_noerr) call handle_err(retval)
C         retval = nf_put_vara_real(ncid, T_varid, start, count, 
C     +        T)
C         if (retval .ne. nf_noerr) call handle_err(retval)

         retval = nf_inq_varid(ncid, 'Q', q_varid)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_put_vara_real(ncid, q_varid, start, count, 
     +        Q)
         if (retval .ne. nf_noerr) call handle_err(retval)

         retval = nf_inq_varid(ncid, 'PSI', psi_varid)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_put_vara_real(ncid, psi_varid, start, count, 
     +        PSI)
         if (retval .ne. nf_noerr) call handle_err(retval)

         retval = nf_inq_varid(ncid, 'U', U_varid)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_put_vara_real(ncid,U_varid, start3, count3, 
     +        U)
         if (retval .ne. nf_noerr) call handle_err(retval)

         retval = nf_inq_varid(ncid, 'QY', QY_varid)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_put_vara_real(ncid,QY_varid, start3, count3, 
     +        QY)
         if (retval .ne. nf_noerr) call handle_err(retval)

C         retval = nf_inq_varid(ncid, 'V', V_varid)
C         if (retval .ne. nf_noerr) call handle_err(retval)
C         retval = nf_put_vara_real(ncid,V_varid,start,count, 
C     +        V)
C         if (retval .ne. nf_noerr) call handle_err(retval)

CCCCC lat-lon only variables
         retval = nf_inq_varid(ncid, 'TS', ts_varid)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_put_vara_real(ncid, ts_varid, start2, count2, 
     +        TS)
         if (retval .ne. nf_noerr) call handle_err(retval)

         retval = nf_inq_varid(ncid, 'TT', tt_varid)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_put_vara_real(ncid, tt_varid, start2, count2, 
     +        TT)
         if (retval .ne. nf_noerr) call handle_err(retval)

CCCCC TIME 
         retval = nf_inq_varid(ncid, 'TIME', tim_varid)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_put_var1_real(ncid, tim_varid, recno,
     +        tim)
         if (retval .ne. nf_noerr) call handle_err(retval)

         retval = nf_inq_varid(ncid, 'GR', gr_varid)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_put_var1_real(ncid, gr_varid, recno,
     +        gr)
         if (retval .ne. nf_noerr) call handle_err(retval)


C     Close the file. This causes netCDF to flush all buffers and make
C     sure your data are really written to disk.
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C non-dimensional time
C      days = tim /86400.d0
C      dys = int(days)
C      hours = (tim - dble(dys)*86400.d0)/3600.d0
C      hrs = int(hours)
C      mins = (tim-dble(dys)*86400.d0-dble(hrs)*3600.d0)/60.d0
C      min = int(mins)

      days = tim/86400.
      dys = int(days)
      hours = (tim - dble(dys)*86400.)/3600.
      hrs = int(hours)
      mins = (tim-dble(dys)*86400.-dble(hrs)*3600.)/60.
      min = int(mins)
      secs = (tim-dble(dys)*86400.-dble(hrs)*3600.-dble(min)*60.)
      sec = int(secs)
   
      write(*,1000) dys,' days '
     +              ,hrs,' hours '
     +              ,min,' minutes '
     +              ,sec,' seconds written to ', FILE_NAME 
C      write(*,2000) ' time '
C     +              ,tim,' written to ', FILE_NAME

 1000 format (I4,A,I2,A,I2,A,I2,A,A)
 2000 format (A,F8.4,A,A)

      end

