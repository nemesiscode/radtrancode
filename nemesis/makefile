#set compiler and flags
FC = $(FCOMP)
FCFLAGS = $(FCFLAGS1) $(STATIC_FLAG)
FRECIPES = $(LIB)/libfrecipes.a


F_SOURCES = \
	coreret.f \
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
	wavesetb.f \
	wavesetc.f \
	mod_scatter.f \
	readkkhead.f \
	readkklblhead.f \
	forwarderr.f \
	forwardavfovX.f \
	forwardavfovMCS.f \
	forwardnogX.f \
	forwardnogXVenus.f \
	forwardnogXrdw.f \
	forwardnoglbl.f \
	forwardnoglblL.f \
	forwardPT.f \
	forwardnogPT.f \
	forwarddisc.f \
	forwardnogL.f \
	forwardnogMCS.f \
	calcnextxn.f \
	setup.f \
	gwritepat.f \
	gwritepatsrom.f \
	extractsrom.f \
	extracttwocloud.f \
	calc_serr.f \
	writeout.f \
	writeoutp.f \
	calcfwderr.f \
	readsurfem.f \
	readraw.f \
	writeraw.f \
        writerawx.f \
	coreretL.f \
	coreretMCS.f \
	forwardavfovL.f \
	gsetrad.f \
	gsetradV.f \
	gsetradL.f \
	gsetradPT.f \
	gsetraddisc.f \
	interpolateMCS.f \
	interpolateMCSnog.f \
	computeFOVA.f \
	computeFOVB.f \
	gsetradMCS.f \
	gwritepatL.f \
	gwritepatV.f \
        gwritepatSO.f \
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
	readcell.f \
	readnextiter.f \
	getradius.f \
	getcloud21.f \
	ackermanmarley.f \
	ackermanmarleyx.f \
	ackermanmarleyx1.f \
	calc_input_files.f \
	tbrownrc.f \
	extractcb.f \
	gwritepatcb.f \
	readmmr.f \
	readrdw.f \
	splang.f \
	splangint.f \
	interpvivien.f \
	modifych4kark.f \
	modifych4irwin.f \
	modifych4sromovsky.f \
	parmentierguillot1.f \
	lineav.f

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
	ar rvU $(LIBRARY) $*.o       
	ranlib $(LIBRARY)
	-rm $*.o


bin :   \
	NemesisX \
	Nemesis \
	NemesisL \
	NemesisPT \
	Nemesisdisc \
	Generatespx \
	GeneratespxL \
	GenerateMCSspx \
	Convertitertomre \
	Testack \
	Testackx \
	Calc_chisq \
	Testbrown \
	Testsrom

NemesisX: nemesis.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/NemesisX nemesis.f $(LIBRARY) \
	$(LIB)/cirsradg.a $(LIB)/cirsrad.a $(LIB)/radtran.a \
	$(LIB)/scatter.a $(LIB)/path.a $(LIB)/matrices.a \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) 

Nemesis: nemesis.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Nemesis nemesis.f $(LIBRARY) \
	$(LIB)/cirsradg.a $(LIB)/cirsrad.a $(LIB)/radtran.a \
	$(LIB)/scatter.a $(LIB)/path.a $(LIB)/matrices.a \
	$(LIB)/spec_data.a $(LIB)/ciatable.a \
	$(LIB)/rtm_util.a $(FRECIPES) 

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

Generatespx: generatespx.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Generatespx generatespx.f $(LIBRARY) \
	$(LIB)/cirsradg.a $(LIB)/cirsrad.a \
	$(LIB)/scatter.a $(LIB)/path.a $(LIB)/radtran.a \
	$(LIB)/matrices.a $(LIB)/spec_data.a $(LIB)/rtm_util.a \
	$(FRECIPES)

GeneratespxL: generatespxL.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/GeneratespxL generatespxL.f $(LIBRARY) \
        $(LIB)/cirsradg.a $(LIB)/cirsrad.a \
        $(LIB)/scatter.a $(LIB)/path.a $(LIB)/radtran.a \
        $(LIB)/matrices.a $(LIB)/spec_data.a $(LIB)/rtm_util.a \
        $(FRECIPES)

PTresults: ptresults.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/PTresults ptresults.f $(LIBRARY) \
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

Convertitertomre: convertitertomre.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Convertitertomre convertitertomre.f \
	$(LIBRARY) $(LIB)/cirsradg.a $(LIB)/cirsrad.a \
	$(LIB)/scatter.a $(LIB)/path.a $(LIB)/radtran.a \
	$(LIB)/matrices.a $(LIB)/spec_data.a $(LIB)/rtm_util.a \
	$(FRECIPES)

Testack: testack.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Testack testack.f \
	$(LIBRARY) $(LIB)/path.a $(LIB)/radtran.a $(LIB)/spec_data.a $(LIB)/rtm_util.a \
	$(FRECIPES)

Testackx: testackx.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Testackx testackx.f \
	$(LIBRARY) $(LIB)/path.a $(LIB)/radtran.a $(LIB)/spec_data.a $(LIB)/rtm_util.a \
	$(FRECIPES)

Testparm: testparm.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Testparm testparm.f \
	$(LIBRARY) $(LIB)/path.a $(LIB)/radtran.a $(LIB)/spec_data.a $(LIB)/rtm_util.a \
	$(FRECIPES)

Teste2: teste2.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Teste2 teste2.f \
	$(LIBRARY) $(LIB)/path.a $(LIB)/radtran.a $(LIB)/spec_data.a $(LIB)/rtm_util.a \
	$(FRECIPES)

Calc_chisq: calc_chisq.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Calc_chisq calc_chisq.f \
	$(LIBRARY) $(LIB)/rtm_util.a $(FRECIPES)

Testbrown: testbrown.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Testbrown testbrown.f \
	$(LIBRARY) $(LIB)/rtm_util.a $(FRECIPES)

Testsrom: testsrom.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Testsrom testsrom.f \
	$(LIBRARY) $(LIB)/path.a $(LIB)/rtm_util.a $(FRECIPES)


