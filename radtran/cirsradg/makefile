############################################################################
#_TITL: MAKEFILE 
#
############################################################################

#set compiler and flags
FC = $(FCOMP)
FCFLAGS = $(FCFLAGS1) $(STATIC_FLAG)
FRECIPES = $(LIB)/libfrecipes.a

F_SOURCES= \
	check_limits.f \
	cirsrtfg_wave.f \
	cirsradg_wave.f \
	close_scat.f \
        get_kg.f \
	map2pro.f \
	map2xvec.f \
	noverlapg.f \
	rankg.f \
	rankk.f \
	sort2g.f \
	subpathg.f \
	calcxmap.f

INCDIR = ../includes
INCLUDES = $(INCDIR)/dbcom.f $(INCDIR)/laycom.f $(INCDIR)/laygrad.f \
	$(INCDIR)/pathcom.f $(INCDIR)/arrdef.f $(INCDIR)/ciacom.f

OBJ = $(F_SOURCES:.f=.o)
LIBRARY = $(LIB)/cirsradg.a
LIBCODE = $(LIBRARY)($(OBJ))


# Target definitions

all: lib

.PHONY: clean
clean:
	-rm -f $(OBJ)/cirsradg.a

.PHONY: touch
touch:
	touch $(F_SOURCES)


lib $(LIBRARY): $(LIBCODE)

$(LIBRARY)(%.o): %.f $(INCLUDES)
	$(FC) $(FCFLAGS) -c $<
	ar rv $(LIBRARY) $*.o
	ranlib $(LIBRARY)
	-rm $*.o


bin :	\
	CIRSdrvg_wave

CIRSdrvg_wave : cirsdrvg_wave.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/CIRSdrvg_wave cirsdrvg_wave.f $(LIBRARY) \
	$(LIB)/cirsrad.a $(LIB)/radtran.a $(LIB)/path.a \
	$(LIB)/scatter.a $(LIB)/rtm_util.a \
	$(LIB)/radtran.a $(LIB)/spec_data.a $(FRECIPES) \
	$(LIB)/matrices.a


