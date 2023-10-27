############################################################################
#_TITL: MAKEFILE 
#
############################################################################

#set compiler and flags
FC = $(FCOMP)
FCFLAGS = $(FCFLAGS1) $(STATIC_FLAG)
FRECIPES = $(LIB)/libfrecipes.a

F_SOURCES = \
	baillynh3.f \
	bboff.f \
	bright.f \
	getmass.f \
	calc_dop_width.f \
	chico2.f \
	burchco2.f \
	hartmannch4.f \
	hartmannch4A.f \
	hartmannch4B.f \
	calc_esum2.f \
	init_lco.f \
	calc_lco.f \
	calc_esum5.f \
	calc_fkdist_wavec.f \
	calc_fkdist_wavec_dp.f \
	calc_lor_width.f \
	calc_mlk1_k.f \
	co2cia.f \
	n2n2cia.f \
	n2h2cia.f \
	calc_mlk2_k.f \
	calc_nlte.f \
	calc_pcoeff.f \
	calc_gvk_k.f \
	nciacon.f \
	ciacon.f \
	co2cia.f \
	e3.f \
	e3interp.f \
	tran_eks.f \
	voigtf_eks.f \
	subline.f \
	linecontrib.f \
	fparacon.f \
	nparacon.f \
	nparacon_all.f \
	ciaread.f \
	ch4cont.f \
	nh3cont.f \
	chico2.f \
	burchco2.f \
	conray.f \
	transet.f \
	roughk.f \
	finefitk.f \
	display.f \
	dpexp.f \
	dslvf.f \
	esumset.f \
	fit_mlband.f \
	gascon.f \
	ngascon.f \
	genrads.f \
	godson_lor.f \
	goody_lor.f \
	gvoico2.f \
	humlic.f \
	intrpk.f \
	fpread.f \
	fpread_all.f \
	kml.f \
	lbl_kcont.f \
	lbl_kdists.f \
	lbl_fknew.f \
	lbl_fknew_dp.f \
	lbl_know.f \
	loadbins.f \
	malkmus_lor.f \
	ml_lacis.f \
	n2n2cia.f \
	n2h2cia.f \
	newconv.f \
	planck.f \
	planckg.f \
	planck_wave.f \
	planckg_wave.f \
	rdlbld.f \
	readdustod.f \
	read_ktable.f \
	read_rodgers.f \
	read_scatter1.f \
	read_xsec.f \
	read_voigt.f \
	slvf.f \
	getnh3.f \
	tau_mg_lor.f \
	tau_goody_voigt1.f \
	tgv3.f \
	tau_eks.f \
	tran_curve.f \
	tran_ml.f \
	voigt_interp.f \
	vvweiss.f \
	tkark.f \
	zgauleg.f \
	getavbins.f \
	getavbinsc.f \
	open_k.f \
	open_kc.f \
	write_k.f \
	reservegas.f \
	readbass1.f \
	interpbass1.f \
	read_ozone_serdyuchenko.f \
	interp_ozone_serdyuchenko.f \
	interp_c2h2_xsections.f \
	interp_c2h6_xsections.f \
	interp_c4h2_xsections.f \
	interp_c3h8_xsections.f \
	interp_ch4_xsections.f \
	interp_co_xsections.f \
	interp_geh4_xsections.f \
	interp_h2o_xsections.f \
	interp_nh3_xsections.f \
	interp_ph3_xsections.f \
	benreuven.f \
	interp_c2h2_xsections.f \
	hmin_ff.f \
	hmin_bf.f


INCDIR = ../includes
INCLUDES = $(INCDIR)/arrdef.f $(INCDIR)/bincom.f $(INCDIR)/bincoms.f \
	$(INCDIR)/dbcom.f $(INCDIR)/ciacom.f \
	$(INCDIR)/lincom.f $(INCDIR)/parcom.f $(INCDIR)/pathcom.f

OBJ = $(F_SOURCES:.f=.o)
LIBRARY = $(LIB)/radtran.a
LIBCODE = $(LIBRARY)($(OBJ))


# Target definitions

all: lib

.PHONY: clean
clean:
	-rm -f $(OBJ)/radtran.a

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
	Aground \
	Ave_table \
	Channel_ave_table \
	Calc_fnktablec \
	Calc_fnktablec_dp \
	Calc_lbltable \
	Aband_ktable \
	Aband_ktablec \
	Conv_spec \
	Convert_table \
	Direct_ktable \
	Direct_ktablec \
	Par_ktable \
	Li_spec \
	Pl_spec \
	Zerotable \
	Read_table \
	Read_lbltable \
	Concat_table \
	Concat_lbl_table \
	Combi_ktable \
	Combi_ktable1 \
	Concat_table_temperature \
	Cut_table \
	Mod_table \
	Mod_table1 \
	Mod_table1a \
	Mod_table1b \
	Mod_table1c \
	Mod_table2 \
	Mod_tableco2 \
	Table_trans 
#	Testsubline
	
#		Radtrans \

	
Testsubline : testsubline.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Testsubline testsubline.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a

Aground : aground.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Aground aground.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a

Ave_table : ave_table.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Ave_table ave_table.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES)  \
	$(LIB)/scatter.a $(LIB)/path.a $(LIB)/matrices.a

Channel_ave_table : channel_ave_table.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Channel_ave_table channel_ave_table.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a $(LIB)/cirsrad.a \
	$(LIB)/scatter.a $(LIB)/path.a $(LIB)/matrices.a $(FRECIPES)

Calc_fnktablec : calc_fnktablec.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Calc_fnktablec calc_fnktablec.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) \
	$(LIB)/scatter.a $(LIB)/path.a $(LIB)/matrices.a

Calc_fnktablec_dp : calc_fnktablec_dp.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Calc_fnktablec_dp calc_fnktablec_dp.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) \
	$(LIB)/scatter.a $(LIB)/path.a $(LIB)/matrices.a
	
Calc_lbltable : calc_lbltable.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Calc_lbltable calc_lbltable.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) \
	$(LIB)/scatter.a $(LIB)/path.a $(LIB)/matrices.a

Zerotable : zerotable.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Zerotable zerotable.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) \
	$(LIB)/scatter.a $(LIB)/path.a $(LIB)/matrices.a

Aband_ktable : aband_ktable.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Aband_ktable aband_ktable.f $(LIBRARY) \
	$(LIB)/rtm_util.a $(LIB)/cirsrad.a $(LIB)/lnemesis.a \
	$(FRECIPES) $(LIB)/matrices.a

Aband_ktablec : aband_ktablec.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Aband_ktablec aband_ktablec.f \
	$(LIBRARY) $(LIB)/rtm_util.a $(LIB)/cirsrad.a $(LIB)/lnemesis.a \
	$(FRECIPES) $(LIB)/matrices.a

Direct_ktable : direct_ktable.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Direct_ktable direct_ktable.f \
	$(LIBRARY) $(LIB)/rtm_util.a \
	$(FRECIPES) $(LIB)/matrices.a

Direct_ktablec : direct_ktablec.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Direct_ktablec direct_ktablec.f \
	$(LIBRARY) $(LIB)/rtm_util.a \
	$(FRECIPES) $(LIB)/matrices.a

Par_ktable : par_ktable.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Par_ktable par_ktable.f \
	$(LIBRARY) $(LIB)/rtm_util.a \
	$(FRECIPES) $(LIB)/matrices.a

Calc_table : calc_table.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Calc_table calc_table.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) \
	$(LIB)/scatter.a $(LIB)/path.a $(LIB)/matrices.a

Conv_spec : conv_spec.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Conv_spec conv_spec.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a

Convert_table : convert_table.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Convert_table convert_table.f \
	$(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) \
	$(LIB)/scatter.a $(LIB)/path.a

Li_spec : li_spec.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Li_spec li_spec.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a

Pl_spec : pl_spec.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Pl_spec pl_spec.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a

Radtrans : radtrans.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Radtrans radtrans.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/scatter.a \
	$(LIB)/path.a $(LIB)/cirsrad.a $(LIB)/rtm_util.a \
	$(LIB)/matrices.a $(FRECIPES) $(LIB)/spec_data.a

Read_table : read_table.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Read_table read_table.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) \
	$(LIB)/scatter.a $(LIB)/path.a $(LIB)/ciatable.a

Read_lbltable : read_lbltable.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Read_lbltable read_lbltable.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) \
	$(LIB)/scatter.a $(LIB)/path.a $(LIB)/ciatable.a

Concat_table : concat_table.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Concat_table concat_table.f \
	$(LIBRARY) $(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) \
	$(LIB)/scatter.a $(LIB)/path.a

Concat_lbl_table : concat_lbl_table.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Concat_lbl_table concat_lbl_table.f \
	$(LIBRARY) $(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) \
	$(LIB)/scatter.a $(LIB)/path.a

Combi_ktable : combi_ktable.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Combi_ktable combi_ktable.f \
	$(LIBRARY) $(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) \
	$(LIB)/scatter.a $(LIB)/path.a

Combi_ktable1 : combi_ktable1.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Combi_ktable1 combi_ktable1.f \
	$(LIBRARY) $(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) \
	$(LIB)/scatter.a $(LIB)/path.a

Concat_table_temperature : concat_table_temperature.f $(LIBRARY)\
	$(LIB)/spec_data.a 
	$(FC) $(FCFLAGS) -o $(BIN)/Concat_table_temperature concat_table_temperature.f \
	$(LIBRARY) $(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) \
	$(LIB)/scatter.a $(LIB)/path.a

Cut_table : cut_table.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Cut_table cut_table.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) \
	$(LIB)/scatter.a $(LIB)/path.a

Mod_table : mod_table.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Mod_table mod_table.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) \
	$(LIB)/scatter.a $(LIB)/path.a

Mod_tableco2 : mod_tableco2.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Mod_tableco2 mod_tableco2.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) \
	$(LIB)/scatter.a $(LIB)/path.a

Mod_table1 : mod_table1.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Mod_table1 mod_table1.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) \
	$(LIB)/scatter.a $(LIB)/path.a

Mod_table1a : mod_table1a.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Mod_table1a mod_table1a.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) \
	$(LIB)/scatter.a $(LIB)/path.a

Mod_table1b : mod_table1b.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Mod_table1b mod_table1b.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) \
	$(LIB)/scatter.a $(LIB)/path.a

Mod_table1c : mod_table1c.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Mod_table1c mod_table1c.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) \
	$(LIB)/scatter.a $(LIB)/path.a

Mod_table2 : mod_table2.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Mod_table2 mod_table2.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) \
	$(LIB)/scatter.a $(LIB)/path.a

Table_trans : table_trans.f $(LIBRARY) $(LIB)/spec_data.a
	$(FC) $(FCFLAGS) -o $(BIN)/Table_trans table_trans.f $(LIBRARY) \
	$(LIB)/spec_data.a $(LIB)/rtm_util.a $(FRECIPES) \
	$(LIB)/scatter.a $(LIB)/path.a

