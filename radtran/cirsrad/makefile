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
	lblconv1.f \
	lblconv2.f \
	cirsrad_wave.f \
	cirsrad_waveS.f \
	cirsrtf.f \
	cirsrtf_wave.f \
	cirsrtf_waveS.f \
	calcsoltau.f \
	close_scat.f \
	findloc.f \
	mradfield.f \
        get_k.f \
        get_klbl.f \
        get_albedo.f \
	get_scatter.f \
	opensol.f \
	init_solar_wave.f \
	get_solar_wave.f \
	get_xsec.f \
	get_xsec_red.f\
	overlap.f \
	rank.f \
	read_khead.f \
	read_klist.f \
	read_klblhead.f \
	read_klbllist.f \
	subpath.f \
	limbscatter.f \
        scatsource.f \
        scatsourcesol.f \
	subcirsrtf_wave.f \
	subcirsrtf_waveS.f \
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
	calc_fine.f \
	readsurfem.f \
	readprfheight.f

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
        CIRSdrv_wavePY \
	Lbldrv_wave

CIRSdrv_wave : cirsdrv_wave.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/CIRSdrv_wave cirsdrv_wave.f $(LIBRARY) \
	$(LIB)/radtran.a $(LIB)/path.a \
	$(LIB)/scatter.a $(LIB)/rtm_util.a $(LIB)/radtran.a \
	$(LIB)/spec_data.a $(LIB)/ciatable.a $(FRECIPES) $(LIB)/matrices.a

CIRSdrv_wavePY : cirsdrv_wavePY.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/CIRSdrv_wavePY cirsdrv_wavePY.f $(LIBRARY) \
        $(LIB)/radtran.a $(LIB)/path.a \
        $(LIB)/scatter.a $(LIB)/rtm_util.a $(LIB)/radtran.a \
        $(LIB)/spec_data.a $(LIB)/lnemesis.a $(LIB)/ciatable.a $(FRECIPES) $(LIB)/matrices.a

Lbldrv_wave : lbldrv_wave.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Lbldrv_wave lbldrv_wave.f $(LIBRARY) \
	$(LIB)/radtran.a $(LIB)/path.a \
	$(LIB)/scatter.a $(LIB)/rtm_util.a $(LIB)/radtran.a \
	$(LIB)/spec_data.a $(FRECIPES) $(LIB)/matrices.a \
	$(LIB)/ciatable.a $(LIB)/rtm_util.a

