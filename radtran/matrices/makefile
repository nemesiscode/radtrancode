############################################################################
#_TITL: MAKEFILE 
#
############################################################################

#set compiler and flags
FC = $(FCOMP)
FCFLAGS = $(FCFLAGS1) $(STATIC_FLAG)
FRECIPES = $(LIB)/libfrecipes.a

F_SOURCES = \
	dadd_mat.f \
	dludcmp.f \
	dlubksb.f \
	dmatinv.f \
	dmult_mat.f \
	dgaussj.f \
	gaussj.f \
	dinvertm.f \
	dcholdc.f \
	dcholsl.f \
	dludcmp1.f \
	gaussjd.f \
	dinvmarq.f \
	dmult_vec.f

INCLUDES = arraycom.f $(INCDIR)/arrdef.f

OBJ= $(F_SOURCES:.f=.o)
LIBRARY = $(LIB)/matrices.a
LIBCODE= $(LIBRARY)($(OBJ))

# Target definitions

all: lib

.PHONY: clean
clean:
	-rm -f $(OBJ)/matrices.a

.PHONY: touch
touch:
	touch $(F_SOURCES)


lib $(LIBRARY): $(LIBCODE)

$(LIBRARY)(%.o): %.f $(INCLUDES)
	$(FC) $(FCFLAGS) -c $<
	ar rv $(LIBRARY) $*.o
	ranlib $(LIBRARY)
	-rm $*.o

