      subroutine ncwrtics(N,M,lons,lats,psi)

      implicit none
      include 'netcdf.inc'

      integer N,M

      integer ncid,retval
      integer lon_dimid, lat_dimid, rec_dimid
      integer lon_varid, lat_varid, rec_varid

      integer dimids(2),start(2),count(2)
      integer psi_varid

      integer ii,jj,tt
      real psi(N,M),pie
      real lons(N),lats(M)

      retval = nf_create('initial.cdf', nf_clobber, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_dim(ncid, 'LON', N, lon_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, 'LAT', M, lat_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, 'lon', NF_REAL, 1, lon_dimid,
     +     lon_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, 'lat', NF_REAL, 1, lat_dimid,
     +     lat_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      dimids(1) = lon_dimid
      dimids(2) = lat_dimid

CCCCC Define the netCDF variable for the data.
CCCCC lat-lon variables
      retval = nf_def_var(ncid, 'INITIAL_PSI', NF_REAL, 2, dimids,
     +     psi_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_var_real(ncid, lon_varid, lons)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, lat_varid, lats)
      if (retval .ne. nf_noerr) call handle_err(retval)

      count(1) = N
      count(2) = M
      start(1) = 1
      start(2) = 1
      retval = nf_inq_varid(ncid, 'INITIAL_PSI', psi_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_real(ncid, psi_varid,
     $			start,count,psi)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
     
      return
      end


