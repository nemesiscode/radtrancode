############################################################################
#_TITL: MAKEFILE 
#
############################################################################

#set compiler and flags
FC = $(FCOMP)
FCFLAGS = $(FCFLAGS1) $(STATIC_FLAG)
FRECIPES = $(LIB)/libfrecipes.a

F_SOURCES = \
	cirsconv.f \
	lblconv.f \
	cirsrad_wave.f \
	cirsrad_waveS.f \
	cirsrad_waveMC.f \
	cirsrtf.f \
	cirsrtf_wave.f \
	cirsrtf_waveS.f \
	cirsrtf_waveMC.f \
	calcsoltau.f \
	close_scat.f \
	findloc.f \
	mradfield.f \
        get_k.f \
        get_albedo.f \
	get_scatter.f \
	opensol.f \
	init_solar_wave.f \
	get_solar_wave.f \
	get_xsec.f \
	overlap.f \
	rank.f \
	read_khead.f \
	read_klist.f \
	subpath.f \
	limbscatter.f \
        scatsource.f \
        scatsourcesol.f \
	subcirsrtf_wave.f \
	subcirsrtf_waveS.f \
	subcirsrtf_waveMC.f \
	dumpflux.f \
	streamflux.f \
	impflux.f \
	dumpfluxsol.f \
	impfluxsol.f \
	openflux.f \
	combineflux.f \
	interpem.f \
	calc_angle.f \
	polar2cartesian.f \
	subview.f \
	init_cont.f \
	calctempparam.f \
	calc_cont.f \
	lblrtf_wave.f \
	lblrad_wave.f \
	calc_fine.f

INCDIR = ../includes
INCLUDES = $(INCDIR)/bincom.f $(INCDIR)/dbcom.f $(INCDIR)/laycom.f \
	$(INCDIR)/pathcom.f $(INCDIR)/arrdef.f $(INCDIR)/ciacom.f

OBJ = $(F_SOURCES:.f=.o)
LIBRARY = $(LIB)/cirsrad.a
LIBCODE = $(LIBRARY)($(OBJ))


# Target definitions

all: lib

.PHONY: clean
clean:
	-rm -f $(OBJ)/cirsrad.a

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
	CIRSdrv_wave \
	Lbldrv_wave

CIRSdrv_wave : cirsdrv_wave.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/CIRSdrv_wave cirsdrv_wave.f $(LIBRARY) \
	$(LIB)/radtran.a $(LIB)/path.a \
	$(LIB)/scatter.a $(LIB)/rtm_util.a $(LIB)/radtran.a \
	$(LIB)/spec_data.a $(FRECIPES) $(LIB)/matrices.a

Lbldrv_wave : lbldrv_wave.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Lbldrv_wave lbldrv_wave.f $(LIBRARY) \
	$(LIB)/radtran.a $(LIB)/path.a \
	$(LIB)/scatter.a $(LIB)/rtm_util.a $(LIB)/radtran.a \
	$(LIB)/spec_data.a $(FRECIPES) $(LIB)/matrices.a \
	$(LIB)/rtm_util.a

