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


part =  ../Apps/pdf.o ../Apps/sfn.o \
        ../theory/altpar15.o ../theory/cteq61pdf.o \
	../theory/xpij2p.o ../theory/xpns2p.o \
	../util/polint.o ../util/alpha_cwz_vmc.o \
	../util/io.o ../util/trmstr.o ../util/linalg.o \
	../nucl/off_shell15.o ../nucl/deuteronSF15.o ../nucl/phi_int2d.o\
	../theory/DIS13.o

all:
	@echo FC = $(FC)
	@echo MATHLIB = $(MATHLIB)

%.o : %.f90
	${FC} -c  $<
	mv $(notdir $@)  util/

PDFerr : PDFerr.o ../Apps/pdf.o $(part)
	$(FC) -g -o ./$@ $(FFLAGS) PDFerr.o ../Apps/pdf.o $(part) $(MATHLIB)

PDFerr18 : PDFerr.o $(part) ../util/alphaS.f90
	$(FC) -g -o ./$@ $(FFLAGS) PDFerr.o $(part) $(MATHLIB) ../util/alphaS.f90

tst_CJpdf: tst_CJpdf.o CJpdf.o
	$(FC) $(FFLAGS) -g -o ./$@ tst_CJpdf.o CJpdf.o

tst_CJsfn: tst_CJsfn.o CJpdf.o
	$(FC) $(FFLAGS) -g -o ./$@ tst_CJsfn.o CJpdf.o


