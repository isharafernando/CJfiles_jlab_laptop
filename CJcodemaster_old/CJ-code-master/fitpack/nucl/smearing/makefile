FC=g77

FFLAGS= -O -fno-automatic -ffixed-line-length-132

PART =  smearfns.o \
	av18.o wjc.o paris.o cdbonn.o \
	../../util/pinterp.o \

smeargrids: smeargrids.f $(PART)
	$(FC) $(FFLAGS) -g -o ./$@ smeargrids.f $(PART)

