############################################################################
#_TITL: MAKEFILE 
#
############################################################################

#set compiler and flags
FC = $(FCOMP)
FCFLAGS = $(FCFLAGS1) $(STATIC_FLAG)
FRECIPES = $(LIB)/libfrecipes.a

F_SOURCES = \
	addgas.f \
	atm.f \
	atmg.f \
	cell.f \
	cloudlayer.f \
        complayer.f \
	cubint.f \
	layer.f \
	layerg.f \
        mewlayer.f \
	muselayer.f \
	newgrav.f \
	newlayer.f \
	uranlayer.f \
	rdmod.f \
	rddmod.f \
	rdamform.f \
	read_grav.f \
	reflatm.f \
	rfpmod.f \
	rfcmod.f \
	sngatm.f \
	verint.f \
	verintg.f \
	calcmolwt.f \
	adjustvmr.f \
	wrlbld.f \
	XhydrostatP.f \
	XhydrostatH.f

INCDIR = ../includes
INCLUDES= $(INCDIR)/pathcom.f $(INCDIR)/laycom.f $(INCDIR)/laygrad.f \
	$(INCDIR)/arrdef.f

OBJ= $(F_SOURCES:.f=.o)
LIBRARY= $(LIB)/path.a
LIBCODE= $(LIBRARY)($(OBJ))

# Target definitions

all: lib

.PHONY: clean
clean:
	-rm -f $(OBJ)/path.a

.PHONY: touch
touch:
	touch $(F_SOURCES)


lib $(LIBRARY): $(LIBCODE)

$(LIBRARY)(%.o): %.f $(INCLUDES)
	$(FC) $(FCFLAGS) -c $<
	ar rv $(LIBRARY) $*.o
	ranlib $(LIBRARY)
	-rm $*.o


bin:     \
	Parah2_profile \
	Dust_profile \
	Nemesisprofile \
	Mod_profile \
	Path \
	Profile \
	Interp_prf \
	Convert_prf \
	Write_path

Dust_profile : dust_profile.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Dust_profile dust_profile.f $(LIBRARY) \
	$(LIB)/radtran.a $(LIB)/spec_data.a $(LIB)/rtm_util.a

Nemesisprofile : nemesisprofile.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Nemesisprofile nemesisprofile.f $(LIBRARY) \
	$(LIB)/rtm_util.a

Parah2_profile : parah2_profile.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Parah2_profile parah2_profile.f $(LIBRARY) \
	$(LIB)/ciatable.a $(LIB)/rtm_util.a

Mod_profile : mod_profile.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Mod_profile mod_profile.f $(LIBRARY) \
	$(LIB)/rtm_util.a

Path : path.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Path path.f $(LIBRARY) \
	$(LIB)/radtran.a $(LIB)/spec_data.a $(LIB)/rtm_util.a

Profile : profile.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Profile profile.f $(LIBRARY) \
	$(LIB)/radtran.a $(LIB)/spec_data.a $(LIB)/rtm_util.a

Interp_prf : interp_prf.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Interp_prf interp_prf.f $(LIBRARY) \
	$(LIB)/rtm_util.a

Convert_prf : convert_prf.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Convert_prf convert_prf.f $(LIBRARY) \
	$(LIB)/rtm_util.a

Write_path : write_path.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Write_path write_path.f \
	$(LIBRARY) $(LIB)/rtm_util.a

