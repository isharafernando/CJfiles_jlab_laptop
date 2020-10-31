#
# File:
#    Makefile
#
# Description:
#    Makefile for compilation of dataDecode program
#
#
DEBUG	?= 1
QUIET	?= 1
#
ifeq ($(QUIET),1)
        Q = @
else
        Q =
endif

CROSS_COMPILE		=
CC			= $(CROSS_COMPILE)gcc -m32
AR                      = ar
RANLIB                  = ranlib
INCS			= -I.
CFLAGS			= -L.
ifeq ($(DEBUG),1)
	CFLAGS		+= -Wall -g
endif
SRC			= $(wildcard *.c)
SRC2			= $(filter-out dataDecode.c,${SRC})
DEPS			= $(SRC2:.c=.d)
PROGS			= $(SRC2:.c=)

all: $(PROGS)

clean distclean:
	@rm -f $(PROGS) *~ $(OBJS) $(DEPS)

%: %.c
	@echo " CC     $@"
	${Q}$(CC) $(CFLAGS) $(INCS) -DDATADECODE=$(@:dec=DataDecode) \
		-include $(@).h -o $@ $< dataDecode.c

%.d: %.c
	@echo " DEP    $@"
	@set -e; rm -f $@; \
	$(CC) -MM -shared $(INCS) -include $(@:.d=.h) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

-include $(DEPS)

.PHONY: all clean distclean
