# Include local system definitions

include ../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)libsubs.a
LIB_MOD=$(LIB_DIR)libmod.a

#
# Entries for library
#
OBJ = $(LIB)(comp_dust_opac.o)\
      $(LIB)(rd_dust_cntrl_file.o)\
      $(LIB)(rd_dust_file.o)\
      $(LIB)(set_ejecta_dust_properties.o)

all : $(LIB)
	rm -f *.o

$(LIB) : $(OBJ)	

DUST = $(LIB_MOD)libmod.a(dust_mod.o)

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f $(DUST)
	$(F90) -c  $(FD) $< 
	ar $(AR_OPTION) $(LIB) $*.o

$(DUST): 
	$(F90) -c  $(FD) dust_mod.f 
	ar $(AR_OPTION) $(LIB_MOD) dust_mod.o

# Specific dependency of source files which uses external 
# fortran moudles. (The order of compilation matters.)


clean :
	rm -f *.o 

clean_lib :
	rm $(LIB)
