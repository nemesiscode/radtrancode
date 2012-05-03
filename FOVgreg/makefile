#set compiler and flags
FC = $(FCOMP)
FCFLAGS = $(FCFLAGS1) $(STATIC_FLAG)
#FCFLAGS = -mcmodel=large -i-dynamic -O -w90 -w95 -FR -cm -w
FCFLAGS = -ffree-form

FRECIPES = $(LIB)/libfrecipes.a

F_SOURCES = \
	FOV_ReadIn.f \
	FOV_Rotation.f \
	verint.f \
	computeFOVgreg.f

OBJ= $(F_SOURCES:.f=.o)
LIBRARY = $(LIB)/lgregfov.a
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
	Test_driver_FOV_nick_pos \
	Outputfov

Test_driver_FOV_nick_pos: test_driver_FOV_nick_pos.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Test_driver_FOV_nick_pos test_driver_FOV_nick_pos.f $(LIBRARY)

Outputfov: outputfov.f $(LIBRARY)
	$(FC) $(FCFLAGS) -o $(BIN)/Outputfov outputfov.f $(LIBRARY)

