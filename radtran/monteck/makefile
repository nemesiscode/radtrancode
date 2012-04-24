#set compiler and flags
FC = $(FCOMP)
FCFLAGS = $(FCFLAGS1) $(STATIC_FLAG)
FRECIPES = $(LIB)/libfrecipes.a


F_SOURCES = \
	arctan.f \
	calcalt.f \
	deflect.f \
	gentabsck.f \
	gentabsck1.f \
	gethg.f \
	hitsun.f \
	henyeymc.f \
	interphg.f \
	modvec.f \
	newpath.f \
	taupath.f \
	phasprob.f \
	readdust.f \
	readprof.f \
	mcphotonck.f \
	yrotate.f

INCDIR= /home/oxpln0/plan/irwin/baldrick/radtran/includes
INCLUDES= 

OBJ = $(F_SOURCES:.f=.o)
LIBRARY = $(LIB)/monteck.a
LIBCODE = $(LIBRARY)($(OBJ))


# target definitions

all: lib


.PHONY: clean
clean:
	-rm -f $(OBJ)/monteck.a

.PHONY: touch
touch:
	touch $(F_SOURCES)


lib $(LIBRARY): $(LIBCODE)

$(LIBRARY)(%.o): %.f $(INCLUDES)
	$(FC) $(FCFLAGS) -c $<
	ar rv $(LIBRARY) $*.o
	ranlib $(LIBRARY)
	-rm $*.o

bin:    \
        Monteck \
	Testnewpath

Monteck : monteck.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Monteck monteck.f \
	$(LIBRARY) $(LIB)/path.a $(LIB)/radtran.a \
	$(LIB)/cirsrad.a $(LIB)/scatter.a \
	$(LIB)/rtm_util.a $(LIB)/radtran.a $(FRECIPES)

Testnewpath : testnewpath.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Testnewpath testnewpath.f \
	$(LIBRARY) $(LIB)/path.a $(LIB)/radtran.a \
	$(LIB)/cirsrad.a $(LIB)/scatter.a \
	$(LIB)/rtm_util.a $(LIB)/radtran.a $(FRECIPES)

