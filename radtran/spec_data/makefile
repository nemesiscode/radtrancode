############################################################################
#_TITL:	MAKEFILE
#
############################################################################

#set compiler and flags
FC = $(FCOMP)
FCFLAGS = $(FCFLAGS1) $(STATIC_FLAG)
FRECIPES = $(LIB)/libfrecipes.a

F_SOURCES = \
	edline.f \
	edset.f \
	edsetnull.f \
	fndwav.f \
	lblband.f \
	liness.f \
	loadbuffer.f \
	partf.f \
	rdgas.f \
	rdiso.f \
	rdkey.f \
	rdkey_band.f \
	rdkey_corrk.f \
	rdline.f \
	printline.f \
	init_partf.f \
	read_bands.f \
	read_del.f \
	read_yam.f \
	write_band.f

INCDIR = ../includes
#INCLUDES = $(INCDIR)/bincom.f $(INCDIR)/bincoms.f $(INCDIR)/dbcom.f \
#	$(INCDIR)/lincom.f $(INCDIR)/pathcom.f $(INCDIR)/arrdef.f

OBJ = $(F_SOURCES:.f=.o)
LIBRARY = $(LIB)/spec_data.a
LIBCODE = $(LIBRARY)($(OBJ))

# Target definitions

all: lib

.PHONY: clean
clean:
	-rm -f $(OBJ)/spec_data.a

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
	Cp_lines \
	Inv_scan \
	Li_lines \
	Makeband \
	Mergeband \
	Modifyband \
	Makedb \
	Makedbloop \
	Merge \
	Merge_multi_ascending \
	Sortlines \
	Pl_lines \
	Pl_elines \
	Scan \
	Select \
	Convertlines \
	Selecttemp \
	Selecttemploop \
	Summary \
	Write_xsec\
	Li_linestemp

Merge : merge.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Merge merge.f $(LIBRARY) \
	$(LIB)/rtm_util.a
	
Merge_multi_ascending : merge_multi_ascending.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Merge_multi_ascending merge_multi_ascending.f $(LIBRARY) \
	$(LIB)/rtm_util.a

Sortlines : sortlines.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Sortlines sortlines.f $(LIBRARY) \
	$(LIB)/rtm_util.a

Addband : addband.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Addband addband.f $(LIBRARY) \
	$(LIB)/rtm_util.a

Cp_lines : cp_lines.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Cp_lines cp_lines.f $(LIBRARY) \
	$(LIB)/rtm_util.a

Li_lines : li_lines.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Li_lines li_lines.f $(LIBRARY) \
	$(LIB)/rtm_util.a

Li_linestemp : li_linestemp.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Li_linestemp li_linestemp.f $(LIBRARY) \
	$(LIB)/path.a $(LIB)/rtm_util.a

Inv_scan : inv_scan.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Inv_scan inv_scan.f $(LIBRARY) \
	$(LIB)/rtm_util.a

Makeband : makeband.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Makeband makeband.f $(LIBRARY) \
	$(LIB)/rtm_util.a

Makedb : makedb.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Makedb makedb.f $(LIBRARY) \
	$(LIB)/rtm_util.a

Makedbloop : makedbloop.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Makedbloop makedbloop.f $(LIBRARY) \
	$(LIB)/rtm_util.a

Scan: scan.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Scan scan.f $(LIBRARY) \
	$(LIB)/rtm_util.a

Mergeband : mergeband.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Mergeband mergeband.f $(LIBRARY) \
	$(LIB)/rtm_util.a

Modifyband : modifyband.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Modifyband modifyband.f $(LIBRARY) \
	$(LIB)/rtm_util.a

Pl_lines : pl_lines.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Pl_lines pl_lines.f $(LIBRARY) \
	$(LIB)/rtm_util.a

Pl_elines : pl_elines.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Pl_elines pl_elines.f $(LIBRARY) \
	$(LIB)/rtm_util.a

Select : select.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Select select.f $(LIBRARY) \
	$(LIB)/rtm_util.a

Convertlines : convertlines.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Convertlines convertlines.f $(LIBRARY) \
	$(LIB)/rtm_util.a

Selecttemp : selecttemp.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Selecttemp selecttemp.f $(LIBRARY) \
	$(LIB)/path.a $(LIB)/rtm_util.a

Selecttemploop : selecttemploop.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Selecttemploop selecttemploop.f $(LIBRARY) \
	$(LIB)/path.a $(LIB)/rtm_util.a

Summary : summary.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Summary summary.f $(LIBRARY) \
	$(LIB)/rtm_util.a

Write_xsec : write_xsec.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Write_xsec write_xsec.f $(LIBRARY) \
	$(LIB)/rtm_util.a
