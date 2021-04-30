#set compiler and flags
FC = $(FCOMP)
FCFLAGS = $(FCFLAGS1) $(STATIC_FLAG)
FRECIPES = $(LIB)/libfrecipes.a

F_SOURCES = \
	adi.f \
	amoeba.f \
	anneal.f \
	avevar.f \
	balanc.f \
	bcucof.f \
	bcuint.f \
	bessi.f \
	bessi0.f \
	bessi1.f \
	bessj.f \
	bessj0.f \
	bessj1.f \
	bessk.f \
	bessk0.f \
	bessk1.f \
	bessy.f \
	bessy0.f \
	bessy1.f \
	beta.f \
	betacf.f \
	betai.f \
	bico.f \
	bksub.f \
	bnldev.f \
	brent.f \
	bsstep.f \
	caldat.f \
	cel.f \
	chder.f \
	chebev.f \
	chebft.f \
	chebpc.f \
	chint.f \
	chsone.f \
	chstwo.f \
	cntab2.f \
	convlv.f \
	correl.f \
	cosft.f \
	covsrt.f \
	crank.f \
	dbrent.f \
	ddpoly.f \
	derf.f \
	des.f \
	desks.f \
	df1dim.f \
	dfpmin.f \
	dgammln.f \
	dgammp.f \
	dgcf.f \
	dgser.f \
	difeq.f \
	eclass.f \
	eclazz.f \
	eigsrt.f \
	el2.f \
	elmhes.f \
	erf.f \
	erfc.f \
	erfcc.f \
	eulsum.f \
	evlmem.f \
	expdev.f \
	f1dim.f \
	factln.f \
	factrl.f \
	fgauss.f \
	fit.f \
	fixrts.f \
	fleg.f \
	flmoon.f \
	four1.f \
	fourn.f \
	fpoly.f \
	frprmn.f \
	ftest.f \
	gamdev.f \
	gammln.f \
	gammp.f \
	gammq.f \
	gasdev.f \
	gauleg.f \
	gaussj.f \
	gcf.f \
	golden.f \
	gser.f \
	hqr.f \
	hunt.f \
	indexx.f \
	irbit1.f \
	irbit2.f \
	jacobi.f \
	julday.f \
	kendl1.f \
	kendl2.f \
	ksone.f \
	kstwo.f \
	laguer.f \
	lfit.f \
	linmin.f \
	locate.f \
	lubksb.f \
	ludcmp.f \
	mdian1.f \
	mdian2.f \
	medfit.f \
	memcof.f \
	midinf.f \
	midpnt.f \
	mmid.f \
	mnbrak.f \
	mnewt.f \
	moment.f \
	mprove.f \
	mrqcof.f \
	mrqmin.f \
	odeint.f \
	pcshft.f \
	pearsn.f \
	piksr2.f \
	piksrt.f \
	pinvs.f \
	plgndr.f \
	poidev.f \
	polcoe.f \
	polcof.f \
	poldiv.f \
	polin2.f \
	polint.f \
	powell.f \
	predic.f \
	probks.f \
	pzextr.f \
	qcksrt.f \
	qgaus.f \
	qromb.f \
	qromo.f \
	qroot.f \
	qsimp.f \
	qtrap.f \
	quad3d.f \
	ran0.f \
	ran1.f \
	ran2.f \
	ran3.f \
	ran4.f \
	rank.f \
	ratint.f \
	realft.f \
	red.f \
	rk4.f \
	rkdumb.f \
	rofunc.f \
	rtbis.f \
	rtflsp.f \
	rtnewt.f \
	rtsafe.f \
	rtsec.f \
	rzextr.f \
	scrsho.f \
	shell.f \
	shoot.f \
	shootf.f \
	simp1.f \
	simp2.f \
	simp3.f \
	simplx.f \
	sinft.f \
	smooft.f \
	sncndn.f \
	solvde.f \
	sor.f \
	sort.f \
	sort2.f \
	sort3.f \
	sparse.f \
	spctrm.f \
	spear.f \
	splie2.f \
	splin2.f \
	spline.f \
	splint.f \
	svbksb.f \
	svdcmp.f \
	svdfit.f \
	svdvar.f \
	toeplz.f \
	tptest.f \
	tqli.f \
	trapzd.f \
	tred2.f \
	tridag.f \
	ttest.f \
	tutest.f \
	twofft.f \
	vander.f \
	zbrac.f \
	zbrak.f \
	zbrent.f \
	zroots.f \
	link.f

OBJ= $(F_SOURCES:.f=.o)
LIBRARY = $(LIB)/libfrecipes.a
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
