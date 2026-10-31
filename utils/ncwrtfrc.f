      subroutine ncwrtfrc(N,M,T,lons,lats,times,eps)

      implicit none
      include 'netcdf.inc'

      integer N,M,T

      integer ncid,retval
      integer lon_dimid, lat_dimid, rec_dimid
      integer lon_varid, lat_varid, rec_varid

      integer dimids(3),start(3),count(3)
      integer eps_varid

      real eps(N,M,T)
      real lons(N),lats(M),times(T)

      retval = nf_create('forcing.cdf', nf_clobber, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_dim(ncid, 'LON', N, lon_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, 'LAT', M, lat_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, 'TIME', T, rec_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, 'LON', NF_REAL, 1, lon_dimid,
     +     lon_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, 'LAT', NF_REAL, 1, lat_dimid,
     +     lat_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, 'TIME', NF_REAL, 1, rec_dimid,
     +     rec_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      dimids(1) = lon_dimid
      dimids(2) = lat_dimid
      dimids(3) = rec_dimid

CCCCC Define the netCDF variable for the data.
CCCCC lat-lon variables
      retval = nf_def_var(ncid, 'EPS', NF_REAL, 3, dimids,
     +     eps_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_var_real(ncid, lat_varid, lats)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, lon_varid, lons)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, rec_varid, times)
      if (retval .ne. nf_noerr) call handle_err(retval)

      count(1) = N
      count(2) = M
      count(3) = T
      start(1) = 1
      start(2) = 1
      start(3) = 1
      retval = nf_inq_varid(ncid, 'EPS', eps_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_real(ncid, eps_varid,
     $			start,count,eps)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      return
      end


