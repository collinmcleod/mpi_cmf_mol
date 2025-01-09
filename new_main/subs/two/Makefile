# Include local system definitions

include ../../../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)libsubs.a
LIB_MOD=$(LIB_DIR)libmod.a

#
# Entries for library
#
OBJ = $(LIB)(steq_ba_two_phot_rate_mpi_v1.o)\
      $(LIB)(get_j_for_two_phot_mpi_v1.o)\
      $(LIB)(two_phot_rate.o)\
      $(LIB)(two_phot_rate_mpi_v1.o)\
      $(LIB)(set_two_phot_v2.o)\
      $(LIB)(set_two_phot_v3.o)\
      $(LIB)(set_two_phot_atm_mpi_v1.o)\
      $(LIB)(two_phot_opac.o)\
      $(LIB)(two_phot_opac_v3.o)\
      $(LIB)(two_phot_opac_mpi_v1.o)\
      $(LIB)(two_phot_var_opac_mpi_v1.o)

all : $(LIB)
	cp -f $(INSTALL_DIR)perm_mod_file jnk.mod
	rm -f *.o
	mv -f *.mod $(MOD_DIR)

$(LIB) : $(OBJ)	

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f $(TWO) $(LIB)(set_kind_module.o)
	$(F90) -c  $(FG) $< 
	ar $(AR_OPTION) $(LIB) $*.o

# Specific dependency of source files which uses external 
# fortran moudles. (The order of compilation matters.)


TWO = $(LIB_MOD)(two_phot_mod.o)
ST = $(LIB_MOD)(steq_data_mod.o)
VAR = $(LIB_MOD)(mod_var_opac_j.o)

$(LIB)(steq_ba_two_phot_rate_v3.o) : steq_ba_two_phot_rate_v3.f $(TWO) $(ST) $(LIB)(set_kind_module.o) 
	$(F90) -c  $(FG) steq_ba_two_phot_rate_v3.f
	ar $(AR_OPTION) $(LIB) steq_ba_two_phot_rate_v3.o

$(LIB)(two_phot_var_opac_mpi_v1.o) : two_phot_var_opac_mpi_v1.f $(TWO) $(VAR) $(LIB)(set_kind_module.o)
	$(F90) -c  $(FG) two_phot_var_opac_mpi_v1.f
	ar $(AR_OPTION) $(LIB) two_phot_var_opac_mpi_v1.o

$(TWO):
	(cd ../../mod_subs ; make)

$(ST):
	(cd ../../mod_subs ; make)

$(VAR):
	(cd ../../mod_subs ; make)

clean :
	rm -f *.o 

clean_lib :
	rm $(LIB)
