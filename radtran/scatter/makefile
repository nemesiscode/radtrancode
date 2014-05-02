############################################################################
#_TITL: MAKEFILE 
#
############################################################################

#set compiler and flags
FC = $(FCOMP)
FCFLAGS = $(FCFLAGS1) $(STATIC_FLAG)
FRECIPES = $(LIB)/libfrecipes.a

F_SOURCES = \
	abend.f \
	add.f \
	addp.f \
	calc_pmat6.f \
	double1.f \
	dmie.f \
	get_hg.f \
	get_theta.f \
	get_input.f \
	get_xsecA.f \
	hansen.f \
	henyey.f \
	fluxint.f \
	fluxintf.f \
	fluxintsolA.f \
	hgaverage.f \
	subhgphas.f \
	locate.f \
	idown.f \
	interp_phase.f \
	iup.f \
	kk_new_sub.f \
	matinv8.f \
	matrices.f \
	miescat.f \
	mrqmin.f \
	mrqcof.f \
	phase1.f \
	phasint2.f \
	rayleigha.f \
	rayleighj.f \
	rayleighv.f \
	read_hg.f \
	read_xsc.f \
	readphlay.f \
	refice.f \
	refindex.f \
	scloud8.f \
	scloud11wave.f \
	scloud12wave.f \
	scloud11flux.f \
	subfithgm.f \
	modmakephase.f

INCDIR = ../includes
#INCLUDES = $(INCDIR)/pathcom.f $(INCDIR)/arrdef.f

OBJ = $(F_SOURCES:.f=.o)
LIBRARY = $(LIB)/scatter.a
LIBCODE = $(LIBRARY)($(OBJ))

# Target definitions

all: lib

.PHONY: clean
clean:
	-rm -f $(OBJ)/scatter.a

.PHONY: touch
touch:
	touch $(F_SOURCES)


lib $(LIBRARY): $(LIBCODE)

$(LIBRARY)(%.o): %.f $(INCLUDES)
	$(FC) $(FCFLAGS) -c $<
	ar rv $(LIBRARY) $*.o
	ranlib $(LIBRARY)
	-rm $*.o

bin :    \
        Makephase \
	Normxsc \
        Readphase

Makephase : makephase.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Makephase makephase.f $(LIBRARY) \
	$(LIB)/rtm_util.a $(FRECIPES)

Normxsc : normxsc.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Normxsc normxsc.f $(LIBRARY) \
	$(LIB)/rtm_util.a

Readphase : readphase.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Readphase readphase.f $(LIBRARY) \
	$(LIB)/rtm_util.a $(FRECIPES)
