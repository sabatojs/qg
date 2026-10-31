      subroutine init_read(N,M,T,psi,eps,FILE_NAME1,FILE_NAME2)

      implicit none
      include 'netcdf.inc'

C     This is the name of the data file we will create.
      character*(*) FILE_NAME1
      character*(*) FILE_NAME2
      integer ncid

      integer NDIMS
      parameter (NDIMS = 3)
      integer N, M, T

C     The start and count arrays will tell the netCDF library where to
C     write our data.
      integer start(NDIMS), count(NDIMS)
      integer start2(NDIMS-1), count2(NDIMS-1)

      integer eps_varid
      integer psi_varid

C     Program variables to hold the data we will write out. We will only
C     need enough space to hold one timestep of data; one record.
      real PSI(N, M)
      real EPS(N, M, T)

C     Error handling.
      integer retval

      count(1) = N
      count(2) = M
      count(3) = T
      start(1) = 1
      start(2) = 1
      start(3) = 1

      count2(1) = N
      count2(2) = M
      start2(1) = 1
      start2(2) = 1

C     open the file. 
      retval = nf_open(FILE_NAME1, nf_write, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

         retval = nf_inq_varid(ncid, 'INITIAL_PSI', psi_varid)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_get_vara_real(ncid, psi_varid, start2, count2, 
     +        PSI)
         if (retval .ne. nf_noerr) call handle_err(retval)

C     Close the file. This causes netCDF to flush all buffers and make
C     sure your data are really written to disk.
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      write(*,*) 'FILE ',FILE_NAME1,' SUCCESSFULLY OPENED AND READ'

C     open the file. 
      retval = nf_open(FILE_NAME2, nf_write, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

         retval = nf_inq_varid(ncid, 'EPS', eps_varid)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_get_vara_real(ncid, eps_varid, start, count, 
     +        EPS)
         if (retval .ne. nf_noerr) call handle_err(retval)

C     Close the file. This causes netCDF to flush all buffers and make
C     sure your data are really written to disk.
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      write(*,*) 'FILE ',FILE_NAME2,' SUCCESSFULLY OPENED AND READ'

      end



