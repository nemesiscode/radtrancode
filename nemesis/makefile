#set compiler and flags
FC = $(FCOMP)
FCFLAGS = $(FCFLAGS1) $(STATIC_FLAG)
FRECIPES = $(LIB)/libfrecipes.a


F_SOURCES = \
	coreret.f \
	coreretMC.f \
	coreretMCMC.f \
	coreretPT.f \
	coreretdisc.f \
	coreret.f \
	write_covariance.f \
	checkfiles.f \
	calc_gain_matrix.f \
	readrefiplan.f \
	readnextspavX.f \
	writenextspavX.f \
	writemvec.f \
	readnextinfo.f \
	readrefhead.f \
	assess.f \
	subprofretg.f \
	calc_phiret.f \
	calc_phiprior.f \
	calc_chi.f \
	writextmp.f \
	readxtmp.f \
	readapriori.f \
	rmodapriori.f \
	modxvecMCMCA.f \
	modxvecMCMC.f \
	wavesetb.f \
	mod_scatter.f \
	readkkhead.f \
	forwarderr.f \
	forwardavfovX.f \
	forwardavfovMCS.f \
	forwardnogX.f \
	forwardnoglbl.f \
	forwardPT.f \
	forwardnogPT.f \
	forwarddisc.f \
	forwardnogMC.f \
	forwardnogL.f \
	forwardnogMCS.f \
	calcnextxn.f \
	setup.f \
	gwritepat.f \
	gwritepatsrom.f \
	extractsrom.f \
	calc_serr.f \
	writeout.f \
	calcfwderr.f \
	readsurfem.f \
	readraw.f \
	writeraw.f \
	coreretL.f \
	coreretMCS.f \
	forwardavfovL.f \
	gsetrad.f \
	gsetradL.f \
	gsetradPT.f \
	gsetraddisc.f \
	interpolateMCS.f \
	interpolateMCSnog.f \
	computeFOVA.f \
	computeFOVB.f \
	gsetradMCS.f \
	gwritepatL.f \
	gwritepatMCS.f \
	gwritepatPT.f \
	gwritepatdisc.f \
	readdrvh.f \
	scankkx.f \
	stripvar.f \
	rankwave.f \
	intradfield.f \
	npvar.f \
	logflag.f \
	check_iteration.f \
	check_profile.f \
	setifix.f \
	readcell.f

INCDIR = ../../radtran/includes
INCLUDES = $(INCDIR)/bincom.f $(INCDIR)/lincom.f $(INCDIR)/pathcom.f \
           $(INCDIR)/parcom.f $(INCDIR)/dbcom.f $(INCDIR)/bincoms.f \
           $(INCDIR)/arrdef.f $(INCDIR)/ciacom.f arraylen.f

OBJ= $(F_SOURCES:.f=.o)
LIBRARY = $(LIB)/lnemesis.a
LIBCODE= $(LIBRARY)($(OBJ))

# target definitions

all: lib

.PHONY: clean
clean:
	-rm -f $(OBJ)

.PHONY: touch
touch:
	touch $(F_SOURCES)


lib $(LIBRARY): $(LIBCODE)

$(LIBRARY)(%.o): %.f $(INCLUDES)
	$(FC) $(FCFLAGS) -c $<
	ar rv $(LIBRARY) $*.o
	ranlib $(LIBRARY)
	-rm $*.o


bin :   \
	Nemesis \
	NemesisL \
	NemesisMCS \
	NemesisMC \
	NemesisMCMC \
	NemesisPT \
	Nemesisdisc \
	Generatespx \
	GenerateMCSspx \
	Radeqm

Nemesis: nemesis.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Nemesis nemesis.f $(LIBRARY) \
	$(LIB)/cirsradg.a $(LIB)/cirsrad.a $(LIB)/radtran.a \
	$(LIB)/scatter.a $(LIB)/path.a $(LIB)/matrices.a \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) 

NemesisMCMC: nemesisMCMC.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/NemesisMCMC nemesisMCMC.f $(LIBRARY) \
	$(LIB)/cirsradg.a $(LIB)/cirsrad.a $(LIB)/radtran.a \
	$(LIB)/scatter.a $(LIB)/path.a $(LIB)/matrices.a \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) 

NemesisPT: nemesisPT.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/NemesisPT nemesisPT.f $(LIBRARY) \
	$(LIB)/cirsradg.a $(LIB)/cirsrad.a \
	$(LIB)/scatter.a  $(LIB)/path.a $(LIB)/radtran.a \
	$(LIB)/matrices.a $(LIB)/spec_data.a $(LIB)/rtm_util.a \
	$(FRECIPES) 

Nemesisdisc: nemesisdisc.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Nemesisdisc nemesisdisc.f $(LIBRARY) \
	$(LIB)/cirsradg.a $(LIB)/cirsrad.a \
	$(LIB)/scatter.a  $(LIB)/path.a $(LIB)/radtran.a \
	$(LIB)/matrices.a $(LIB)/spec_data.a \
	$(LIB)/rtm_util.a $(FRECIPES)

NemesisL: nemesisL.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/NemesisL nemesisL.f $(LIBRARY) \
	$(LIB)/cirsradg.a $(LIB)/cirsrad.a \
	$(LIB)/scatter.a $(LIB)/path.a $(LIB)/matrices.a \
	$(LIB)/radtran.a $(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) 

NemesisMCS: nemesisMCS.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/NemesisMCS nemesisMCS.f $(LIBRARY) \
	$(LIB)/cirsradg.a $(LIB)/cirsrad.a \
	$(LIB)/scatter.a $(LIB)/path.a $(LIB)/matrices.a \
	$(LIB)/radtran.a $(LIB)/spec_data.a $(LIB)/rtm_util.a \
	$(LIB)/lgregfov.a $(FRECIPES) 

NemesisMC: nemesisMC.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/NemesisMC nemesisMC.f $(LIBRARY) \
	$(LIB)/cirsradg.a $(LIB)/cirsrad.a \
	$(LIB)/monteck.a \
	$(LIB)/scatter.a $(LIB)/path.a $(LIB)/radtran.a \
	$(LIB)/matrices.a $(LIB)/spec_data.a $(LIB)/rtm_util.a \
	$(FRECIPES) 

Generatespx: generatespx.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Generatespx generatespx.f $(LIBRARY) \
	$(LIB)/cirsradg.a $(LIB)/cirsrad.a \
	$(LIB)/scatter.a $(LIB)/path.a $(LIB)/radtran.a \
	$(LIB)/matrices.a $(LIB)/spec_data.a $(LIB)/rtm_util.a \
	$(FRECIPES)

GenerateMCSspx: generateMCSspx.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/GenerateMCSspx generateMCSspx.f \
	$(LIBRARY) \
	$(LIB)/cirsradg.a $(LIB)/cirsrad.a \
	$(LIB)/scatter.a $(LIB)/path.a $(LIB)/radtran.a \
	$(LIB)/matrices.a $(LIB)/spec_data.a $(LIB)/rtm_util.a \
	$(LIB)/lgregfov.a $(FRECIPES)

Radeqm: radeqm.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Radeqm radeqm.f $(LIBRARY) \
	$(LIB)/cirsrad.a \
	$(LIB)/scatter.a $(LIB)/path.a $(LIB)/radtran.a \
	$(LIB)/matrices.a $(LIB)/spec_data.a $(LIB)/rtm_util.a \
	$(FRECIPES)

