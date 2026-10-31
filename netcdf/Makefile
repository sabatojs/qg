NCSRCS = ncinit.f ncwritehist.f ncinit_lin.f ncwritehist_lin.f ncreadhist.f ncwritespec.f

all: $(NCSRCS)
	$(FC) $(NCINCLUDES) -c $(NCSRCS) 

clean:
	rm -rf *.o
