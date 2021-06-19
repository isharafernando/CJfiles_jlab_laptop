# Selects the right gfortran compiler call
ifneq (,$(findstring Bilbo,$(shell uname -n)))
  FC = gfortran-4.7	# for Alberto's odl laptop
else
  FC = gfortran		# for everybody else
endif

FFLAGS= -O -fno-automatic -ffixed-line-length-132

# Selects the right mathlib library
ifneq (,$(findstring _64,$(shell uname -m)))
  MATHLIB = ../util/libmathlib64.a	# for for 64 bits machines
else
  MATHLIB = ../util/libmathlib.a	# for for 32 bits machines
endif


all:
	@echo FC = $(FC)
	@echo MATHLIB = $(MATHLIB)

%.o : %.f90
	${FC} -c  $<
	mv $(notdir $@)  util/


### F2 stand-alone calculations

part =  pdf.o sfn.o\
	../theory/altpar15.o ../theory/xpij2p.o \
	../theory/xpns2p.o \
	../util/alpha_cwz_vmc.o $(MATHLIB)\
	../util/io.o ../util/trmstr.o \

F2c =   ../theory/DIS13.o ../theory/cteq61pdf.o $(part) \
	../nucl/deuteronSF15.o ../nucl/phi_int2d.o ../nucl/phi_int3d.o \
	../nucl/off_shell15.o \
	../util/linalg.f ../util/alphaS.f90 

F2comp : F2comp.o $(F2c) 
	$(FC) -g -o ./$@ F2comp.o $(F2c) 


### Truncated moment calcs

truncmom: truncmom.o $(F2c)
	$(FC) -g -o ./$@ truncmom.o $(F2c)



#............ OLDER STUFF (was in parent folder, check paths) ............

############   Applications compilation and release  #############


#part =  Apps/pdf.o \
#	theory/altpar15.o theory/xpij2p.o \
#	theory/xpns2p.o \
#	util/polint.o util/alpha_cwz.o $(MATHLIB)\
#	util/io.o util/trmstr.o \

#F2c =   theory/DIS13.o $(part) \
#	nucl/deuteronSF.new.o nucl/phi_int2d.o nucl/phi_int3d.o \
#	util/linalg.f	

F2c.r = Apps/F2comp.f Apps/pdf.f makefile_F2comp README.F2comp \
	theory/DIS13.f theory/altpar15.f theory/xpij2p.f theory/xpns2p.f \
	util/polint.f util/alpha_cwz.f $(MATHLIB) util/io.f \
	util/trmstr.f \
	nucl/deuteronSF.f nucl/phi_int2d.f nucl/phi_int3d.f \

FILES_F2c = Apps/CTEQ6X_wba.par Apps/CTEQ6X_wba_mst.par Apps/CTEQ6X_dmc.par  \
	Apps/CTEQ6X_free.par \
	Apps/CTEQ6X_wba.gf Apps/CTEQ6X_wba_mst.gf Apps/CTEQ6X_dmc.gf  \
	Apps/CTEQ6X_free.gf \
	nucl/phi.* 

#TM = 	$(F2c) \
#	util/inttochar.o


#F2comp : Apps/F2comp.o $(F2c) 
#	$(FC) -g -o ./$@ Apps/F2comp.o $(F2c) 

FFcomp : Apps/FFcomp.o $(F2c) 
	$(FC) -g -o ./$@ Apps/FFcomp.o $(F2c) 

F2-BONUS : Apps/F2-BONUS.o $(F2c) 
	$(FC) -g -o ./$@ Apps/F2-BONUS.o $(F2c) 

partlumi : Apps/partlumi.o $(part) 
	$(FC) -g -o ./$@ Apps/partlumi.o $(part) 


F2comp_RELEASE:	
	zip F2comp_RELEASE.zip $(F2c.r) $(FILES_F2c) Apps/F2comp.out_test


#	Test of F2comp [AA, 17 May 10]

test_F2c:	
	make F2comp
	mv F2comp tests/F2comp/
	cd tests/F2comp/ && ./F2comp > F2comp.out  \
			 && kdiff3 F2comp.out F2comp.out_test
