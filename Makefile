# just in case i screw up the makefile this is what i did to compile at the command line
#f77 shallow_water.f  gauss_legendre.f spectral_transform.f Gnm_R_S.f
#f77 -c ncwrapper.f ncinit.f ncwritehist.f -I/home/jssabato/netcdf-3.6.1/include
#f77 -o nctest.x ncinit.o ncwrapper.o ncwritehist.o -L/home/jssabato/netcdf-3.6.1/lib/lib -lnetcdf

EXE_NAME=qg_demo
FC=gfortran
ALL_OBJS = `find . -name "*.o"`
HOMEDIR = /home/jssabato
QGDIR = $(HOMEDIR)/FORTRAN-CODES/qg
NCINCLUDES=-I/usr/lib/netcdf_64b_gnu4.3/include
NETCDF=-L/usr/lib/netcdf_64b_gnu4.3/lib -lnetcdff -lnetcdf
FFLAGS = 

export HOMEDIR
export NCINCLUDES
export FFLAGS
export FC

OBJS   = qg.o

all: $(OBJS) net_cdf utilities
	$(FC) $(FFLAGS) -o $(EXE_NAME) $(ALL_OBJS) $(NETCDF)

utilities:
	cd utils; $(MAKE) all

net_cdf:
	cd netcdf; $(MAKE) all

clean: 
	if [ -f $(EXE_NAME) ]; then rm $(EXE_NAME); fi
	rm $(ALL_OBJS)
