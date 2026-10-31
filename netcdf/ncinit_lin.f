C     This is part of the netCDF package.
C     Copyright 2006 University Corporation for Atmospheric Research/Unidata.
C     See COPYRIGHT file for conditions of use.

      subroutine ncinit_lin(FILE_NAME,NLEVS,NLATS,NLONS,levs,lats,lons,
     $			TH,N2,H,BETA,DEL,f,g,U0,Lx,Ly,gr,kx)
      implicit none
      include 'netcdf.inc'

C     This is the name of the data file we will create.
      character*(*) FILE_NAME
      integer ncid,i,j

C     We are writing 3D data, a M x N lat-lon grid, with
C     unlimited timesteps of data.
      integer NDIMS, NRECS
      parameter (NDIMS = 4)
      integer NLATS, NLONS, NLEVS
      character*(*) LAT_NAME, LON_NAME, LEV_NAME, REC_NAME
      parameter (LAT_NAME = 'LAT', LON_NAME = 'LON',LEV_NAME='LEV')
      parameter (REC_NAME = 'TIME')
      integer lon_dimid, lev_dimid,lat_dimid, rec_dimid

C     The start and count arrays will tell the netCDF library where to
C     write our data.
      integer start(NDIMS), count(NDIMS)
      real TH(NLEVS),N2,H,BETA,f,g,U0,Lx,Ly,DEL,gr,kx

C constants
      integer N2_varid,H_varid,BETA_varid,f_varid,g_varid,
     +        U0_varid,Lx_varid,Ly_varid,del_varid,gr_varid,
     +	      kx_varid

C     These program variables hold the latitudes and longitudes.
      real lats(NLATS), lons(NLONS), levs(NLEVS)
      integer lon_varid, lev_varid, lat_varid, rec_varid

C     create netCDF variables

CCCCC lat-lon-lev variables
      character*(*) U_NAME
      parameter (U_NAME='U')
      integer u_varid
      character*(*) QY_NAME
      parameter (QY_NAME='QY')
      integer QY_varid
      character*(*) TH_NAME
      parameter (TH_NAME='T0')
      integer TH_varid
      character*(*) q_NAME
      parameter (q_NAME='Q')
      integer q_varid
      character*(*) PSI_NAME
      parameter (PSI_NAME='PSI')
      integer PSI_varid

CCCCC lat-lon only variables
      character*(*) TT_NAME
      parameter (TT_NAME='TT')
      integer tt_varid
      character*(*) TS_NAME
      parameter (TS_NAME='TS')
      integer ts_varid

      integer dimids(NDIMS)
      integer dimids2(NDIMS-1)
      integer dimids3(NDIMS-1)

C     We recommend that each variable carry a "units" attribute.
      character*(*) UNITS
      parameter (UNITS = 'units')
      character*(*) T_UNITS, Q_UNITS, U_UNITS, PSI_UNITS,
     +		N2_UNITS,BETA_UNITS,f_UNITS,H_UNITS,g_UNITS,
     +		REC_UNITS,LAT_UNITS,LON_UNITS,LEV_UNITS,
     +		DEL_UNITS,QY_UNITS,gr_UNITS,kx_UNITS
      parameter (T_UNITS = 'K')
      parameter (PSI_UNITS = 'm2/s')
      parameter (Q_UNITS = '1/s')
      parameter (QY_UNITS = '1/ms')
      parameter (U_UNITS = 'm/s')
      parameter (LAT_UNITS = 'm')
      parameter (LON_UNITS = 'm')
      parameter (LEV_UNITS = 'm')
      parameter (REC_UNITS = 's')
      parameter (g_UNITS = 'm/s2')
      parameter (H_UNITS = 'm')
      parameter (f_UNITS = '1/s')
      parameter (N2_UNITS = '1/s2')
      parameter (BETA_UNITS = '1/ms')
      parameter (gr_UNITS = '1/s')
      parameter (kx_UNITS = 'ND')
      parameter (DEL_UNITS = '1/m2s')

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
      retval = nf_def_dim(ncid, LAT_NAME, NLATS, lat_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, LON_NAME, NLONS, lon_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, LEV_NAME, NLEVS, lev_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, REC_NAME, NF_UNLIMITED, rec_dimid)
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

      retval = nf_def_var(ncid, LON_NAME, NF_REAL, 1, lon_dimid, 
     +     lon_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, LEV_NAME, NF_REAL, 1, lev_dimid, 
     +     lev_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, TH_NAME, NF_REAL, 1, lev_dimid, 
     +     TH_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, REC_NAME, NF_REAL, 1, rec_dimid, 
     +     rec_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C constants
      retval = nf_def_var(ncid, 'g', NF_REAL, 0, 0, 
     +     g_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, 'N2', NF_REAL, 0, 0, 
     +     N2_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, 'H', NF_REAL, 0, 0, 
     +     H_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, 'f', NF_REAL, 0, 0, 
     +     f_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, 'BETA', NF_REAL, 0, 0, 
     +     BETA_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, 'GR', NF_REAL, 1, rec_dimid, 
     +     gr_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C      retval = nf_def_var(ncid, 'GR', NF_REAL, 0, 0, 
C     +     gr_varid)
C      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, 'KX', NF_REAL, 0, 0, 
     +     kx_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, 'DELTA', NF_REAL, 0, 0, 
     +     DEL_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, 'U0', NF_REAL, 0, 0, 
     +     U0_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, 'Lx', NF_REAL, 0, 0, 
     +     Lx_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, 'Ly', NF_REAL, 0, 0, 
     +     Ly_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

CCCCCCCC Assign units attributes to coordinate variables.

      retval = nf_put_att_text(ncid, lat_varid, UNITS, len(LAT_UNITS), 
     +     LAT_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, lon_varid, UNITS, len(LON_UNITS), 
     +     LON_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, lev_varid, UNITS, len(LEV_UNITS), 
     +     LEV_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, rec_varid, UNITS, len(REC_UNITS), 
     +     REC_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

C constants
      retval = nf_put_att_text(ncid, g_varid, UNITS, len(G_UNITS), 
     +     G_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, H_varid, UNITS, len(H_UNITS), 
     +     H_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, BETA_varid, UNITS,len(BETA_UNITS), 
     +     BETA_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, gr_varid, UNITS,len(gr_UNITS), 
     +     gr_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, kx_varid, UNITS,len(kx_UNITS), 
     +     kx_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, DEL_varid, UNITS,len(DEL_UNITS), 
     +     DEL_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, N2_varid, UNITS, len(N2_UNITS), 
     +     N2_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, f_varid, UNITS, len(f_UNITS), 
     +     f_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, Lx_varid, UNITS, len(LAT_UNITS), 
     +     LAT_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, Ly_varid, UNITS, len(LAT_UNITS), 
     +     LAT_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, U0_varid, UNITS, len(U_UNITS), 
     +     U_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, TH_varid, UNITS, len(T_UNITS), 
     +     T_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     The dimids array is used to pass the dimids of the dimensions of
C     the netCDF variables. Both of the netCDF variables we are creating
C     share the same four dimensions. In Fortran, the unlimited
C     dimension must come last on the list of dimids.
      dimids(1) = lon_dimid
      dimids(2) = lat_dimid
      dimids(3) = lev_dimid
      dimids(4) = rec_dimid

      dimids2(1)=lon_dimid
      dimids2(2)=lat_dimid
      dimids2(3)=rec_dimid

      dimids3(1)=lat_dimid
      dimids3(2)=lev_dimid
      dimids3(3)=rec_dimid

CCCCC Define the netCDF variables for the data.
CCCCC lat-lon-lev variables

C      retval = nf_def_var(ncid, T_NAME, NF_REAL, NDIMS, dimids, 
C     +     t_varid)
C      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, Q_NAME, NF_REAL, NDIMS, dimids, 
     +     q_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, PSI_NAME, NF_REAL, NDIMS, dimids, 
     +     PSI_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, U_NAME, NF_REAL, NDIMS-1, dimids3, 
     +     U_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, QY_NAME, NF_REAL, NDIMS-1, dimids3, 
     +     QY_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

CCCCC lat-lon only variables
      retval = nf_def_var(ncid, TS_NAME, NF_REAL, NDIMS-1, dimids2, 
     +     ts_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, TT_NAME, NF_REAL, NDIMS-1, dimids2, 
     +     tt_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     Assign units attributes to the netCDF variables.
CCCC lat-lon-lev variables
C      retval = nf_put_att_text(ncid, t_varid, UNITS, len(T_UNITS), 
C     +     T_UNITS)
C      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, q_varid, UNITS,len(Q_UNITS), 
     +     Q_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, PSI_varid, UNITS,len(PSI_UNITS), 
     +     PSI_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, U_varid, UNITS, len(U_UNITS), 
     +     U_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, QY_varid, UNITS,len(QY_UNITS), 
     +     QY_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

CCCC lat-lon only variables
      retval = nf_put_att_text(ncid, ts_varid, UNITS, len(T_UNITS), 
     +     T_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, tt_varid, UNITS, len(T_UNITS), 
     +     T_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     End define mode.
      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     Write the coordinate variable data. This will put the latitudes
C     and longitudes of our data grid into the netCDF file.
      retval = nf_put_var_real(ncid, lat_varid, lats)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, lon_varid, lons)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, lev_varid, levs)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_var_real(ncid, TH_varid, TH)
      if (retval .ne. nf_noerr) call handle_err(retval)


C constants
      retval = nf_put_var_real(ncid, g_varid, g)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, H_varid, H)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, BETA_varid, BETA)
      if (retval .ne. nf_noerr) call handle_err(retval)
C      retval = nf_put_var_real(ncid, gr_varid, gr)
C      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, gr_varid, gr)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, kx_varid, kx)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, DEL_varid, DEL)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, N2_varid, N2)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, f_varid, f)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, U0_varid, U0)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, Lx_varid, Lx)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, Ly_varid, Ly)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
   
      write(*,*)
      print *,'*** SUCCESS CREATING OUTPUT FILE ', FILE_NAME, '!'
      write(*,*)

      end

