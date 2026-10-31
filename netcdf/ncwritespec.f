C     This is part of the netCDF package.
C     Copyright 2006 University Corporation for Atmospheric Research/Unidata.
C     See COPYRIGHT file for conditions of use.

      subroutine ncwritespec(UVSPEC,dUVdySPEC,VQSPEC,k,w,LAT,N,M
     $		,Tf,FILE_NAME)
      implicit none
      include 'netcdf.inc'

C     This is the name of the data file we will create.
      character*(*) FILE_NAME
      integer ncid,i,j,t

C     We are writing 3D data, a M x N lat-lon grid, with
C     unlimited timesteps of data.
      integer NDIMS, NRECS
      parameter (NDIMS = 3)
      integer N, M, Tf
      character*(*) LAT_NAME, k_NAME, w_NAME
      parameter (LAT_NAME = 'LAT', k_NAME = 'WAVE',w_NAME='FREQ')
      integer lat_dimid, k_dimid, w_dimid

C     The start and count arrays will tell the netCDF library where to
C     write our data.
      integer start(NDIMS), count(NDIMS)
      real UVSPEC(N/2,Tf,M),dUVdySPEC(N/2,Tf,M),VQSPEC(N/2,Tf,M)

C     These program variables hold the latitudes and longitudes.
      real LAT(M), k(N/2), w(Tf)
      integer lat_varid, k_varid, w_varid

C     create netCDF variables

CCCCC lat-k-w variables
      character*(*) UV_NAME
      parameter (UV_NAME='UVSPEC')
      integer uv_varid
      character*(*) dUVdy_NAME
      parameter (dUVdy_NAME='dUVdySPEC')
      integer duvdy_varid
      character*(*) VQ_NAME
      parameter (VQ_NAME='VQSPEC')
      integer VQ_varid

      integer dimids(NDIMS)

C     Program variables to hold the data we will write out. We will only
C     need enough space to hold one timestep of data; one record.

C     Error handling.
      integer retval

C     Create the file. 
      retval = nf_create(FILE_NAME, nf_clobber, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     Define the dimensions. The record dimension is defined to have
C     unlimited length - it can grow as needed. In this example it is
C     the time dimension.
      retval = nf_def_dim(ncid, LAT_NAME, M, lat_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, k_NAME, N/2, k_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, w_NAME, Tf, w_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     Define the coordinate variables. We will only define coordinate
C     variables for lat, lon and lev.  Ordinarily we would need to provide
C     an array of dimension IDs for each variable's dimensions, but
C     since coordinate variables only have one dimension, we can
C     simply provide the address of that dimension ID (lat_dimid) and
C     similarly for (lon_dimid) and (lev_dimid).

      retval = nf_def_var(ncid, LAT_NAME, NF_REAL, 1, lat_dimid, 
     +     lat_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, k_NAME, NF_REAL, 1, k_dimid, 
     +     k_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, w_NAME, NF_REAL, 1, w_dimid, 
     +     w_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     The dimids array is used to pass the dimids of the dimensions of
C     the netCDF variables. Both of the netCDF variables we are creating
C     share the same four dimensions. In Fortran, the unlimited
C     dimension must come last on the list of dimids.
      dimids(1) = k_dimid
      dimids(2) = w_dimid
      dimids(3) = lat_dimid

CCCCC Define the netCDF variables for the data.
CCCCC k-w-lat variables
      retval = nf_def_var(ncid, VQ_NAME, NF_REAL, NDIMS, dimids, 
     +     vq_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, UV_NAME, NF_REAL, NDIMS, dimids, 
     +     UV_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, dUVdy_NAME, NF_REAL, NDIMS, dimids, 
     +     dUVdy_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     End define mode.
      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     Write the coordinate variable data. This will put the latitudes
C     and longitudes of our data grid into the netCDF file.
      retval = nf_put_var_real(ncid, lat_varid, LAT)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, k_varid, k)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, w_varid, w)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
   
      write(*,*)
      print *,'*** SUCCESS CREATING OUTPUT FILE ', FILE_NAME, '!'
      write(*,*)

C     open the file. 
      retval = nf_open(FILE_NAME, nf_write, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     These settings tell netcdf to write one timestep of data. (The
C     setting of start(3) from the time step input tells netCDF which
C     timestep to write.)
      count(1) = N/2
      count(2) = Tf
      count(3) = M
      start(1) = 1
      start(2) = 1
      start(3) = 1

CCCCC k-w-lat variables
         retval = nf_inq_varid(ncid, 'VQSPEC', vq_varid)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_put_vara_real(ncid, vq_varid, start, count, 
     +        VQSPEC)
         if (retval .ne. nf_noerr) call handle_err(retval)

         retval = nf_inq_varid(ncid, 'UVSPEC', uv_varid)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_put_vara_real(ncid, uv_varid, start, count, 
     +        UVSPEC)
         if (retval .ne. nf_noerr) call handle_err(retval)

         retval = nf_inq_varid(ncid, 'dUVdySPEC', duvdy_varid)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_put_vara_real(ncid, duvdy_varid, start, count, 
     +        dUVdySPEC)
         if (retval .ne. nf_noerr) call handle_err(retval)

C     Close the file. This causes netCDF to flush all buffers and make
C     sure your data are really written to disk.
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)



      end

