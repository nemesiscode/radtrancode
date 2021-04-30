############################################################################
#_TITL: MAKEFILE 
#
############################################################################

#set compiler and flags
FC = $(FCOMP)
FCFLAGS = $(FCFLAGS1) $(STATIC_FLAG)
FRECIPES = $(LIB)/libfrecipes.a

F_SOURCES = \
        askyn.f \
        bit.f \
        clp.f \
	datarchive.f \
        ebcasc.f \
        simp_int.f \
        fcat.f \
        file.f \
        readflags.f \
        hamming.f \
        hanning.f \
        gasdev1.f \
        getquad.f \
	interp.f \
        locase.f \
        prompt.f \
        ran11.f \
        readi.f \
        readr.f \
        remsp.f \
        upcase.f \
        wtext.f \
        cspline.f \
        csplint.f \
        isys.f \
	bilinear.f

OBJ = $(F_SOURCES:.f=.o)
LIBRARY = $(LIB)/rtm_util.a
LIBCODE = $(LIBRARY)($(OBJ))


INCDIR = ../includes
INCLUDES = $(INCDIR)/arrdef.f

# Target definitions

all: $(LIBRARY)

.PHONY: clean
clean:
	-rm -f $(OBJ)/rtm_util.a

.PHONY: touch
touch:
	touch $(F_SOURCES)


lib $(LIBRARY): $(LIBCODE)

$(LIBRARY)(%.o): %.f $(INCLUDES)
	$(FC) $(FCFLAGS) -c $<
	ar rvU $(LIBRARY) $*.o
	ranlib $(LIBRARY)
	-rm $*.o

bin:    \
