# Include local system definitions

include ../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)libsubs.a
LIB_MOD=$(LIB_DIR)libmod.a


# Entries for library

OBJ = $(LIB)(write_recom_mpi_v1.o)\
      $(LIB)(set_idisp_isend.o)\
      $(LIB)(wr_cool_mpi_v1.o)\
      $(LIB)(wr2d_mpi_v1.o)\
      $(LIB)(wr2d_gath_mpi_v1.o)\
      $(LIB)(wrrecomchk_mpi_v1.o)

all : $(LIB)
	rm -f *.o
$(LIB) : $(OBJ)	

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f $(TWO)
	$(F90) -c  $(FG) $< 
	ar $(AR_OPTION) $(LIB) $*.o

# Specific dependency of source files which uses external 
# fortran moudles. (The order of compilation matters.)


CMF_MOD = $(LIB_MOD)(mod_cmfgen.o)

ST = $(LIB_MOD)(steq_data_mod.o)

$(LIB)(steq_ba_two_phot_rate_v3.o) : steq_ba_two_phot_rate_v3.f $(TWO) 
	$(F90) -c  $(FG) write_recom_mpi_v1.f
	ar $(AR_OPTION) $(LIB) write_recom_mpi_v1.o

$(CMF_MOD):
	(cd ../new_main/mod_subs ; make)

$(ST):
	(cd ../../mod_subs ; make)

clean :
	rm -f *.o 

clean_lib :
	rm $(LIB)
