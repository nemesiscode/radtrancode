############################################################################
#_TITL: MAKEFILE 
#
#_SRCE:	/radtran/ciatable/makefile
#
#_HIST:	3oct02	PDP	Created this header. Modified the file to allow
#			compilation with an Intel FORTRAN compiler under
#			Linux and/or optimation if using a PentiumIII
#			chip. Added n2n2_s.f to the list of F_SOURCES as
#			it was previously missing.
############################################################################

#set compiler and flags
FC = $(FCOMP)
FCFLAGS = $(FCFLAGS1) $(STATIC_FLAG)
FRECIPES = $(LIB)/libfrecipes.a

F_SOURCES = \
      nh2h2.f \
      nh2he.f \
      npartsum.f \
      calcpara.f \
      calcfracpara.f \
      ch4ch4_s.f \
      h2ch4_s.f \
      h2ch4_sub.f \
      h2h2_v0s.f \
      h2h2_v1s.f \
      h2h2_v2s.f \
      h2he_v0s.f \
      h2he_v1s.f \
      n2ch4_s.f \
      n2h2_s.f \
      n2h2_sub.f \
      n2n2_s.f \
      pindex.f

INCDIR = ../includes
INCLUDES = $(INCDIR)/ciacom.f


OBJ = $(F_SOURCES:.f=.o)
LIBRARY = $(LIB)/ciatable.a
LIBCODE = $(LIBRARY)($(OBJ))


# Target definitions

all: lib

.PHONY: clean
clean:
	-rm -f $(OBJ)/ciatable.a

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
	Makeciatable \
	Makeisotable \
	Makefptable \
	Makefptable_allcia \
	Testfracpara \
	Testh2h2v2s

Makeciatable: makeciatable.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Makeciatable makeciatable.f $(LIBRARY)

Makeisotable: makeisotable.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Makeisotable makeisotable.f $(LIBRARY)

Makefptable: makefptable.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Makefptable makefptable.f $(LIBRARY)

Makefptable_allcia: makefptable_allcia.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Makefptable_allcia makefptable_allcia.f $(LIBRARY)


Testfracpara: testfracpara.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Testfracpara testfracpara.f $(LIBRARY)

Testh2h2v2s: testh2h2v2s.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Testh2h2v2s testh2h2v2s.f $(LIBRARY)

Testcalcpara: testcalcpara.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Testcalcpara testcalcpara.f $(LIBRARY)

