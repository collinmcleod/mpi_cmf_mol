# TEST MAKEFILE


# Set parameters for arechiv .o files into libraries.
# On UBUNTU systems msy nrrf Uruv

AR_OPTION = ruv

HOST :=$(shell hostname)

#Include local system dependent definitions

include Makefile_definitions


# We will access the makefile in each local directory tp create the
# libraries and executables.

all : info d_tools d_blas d_dust d_lpack d_nm_subs d_stark d_unix d_subs d_plane d_newsubs \
         d_mpi d_new_main 
     

# We now MAKE the required libraries and executables.

info:
	@echo " "
	@echo " Host         = "$(HOST)
	@echo " Compile opt. = "$(COMPILE_OPT)
	@echo " Install dir. = "$(INSTALL_DIR)
	@echo " (F90)        = "$(F90)
	@echo " (FG)         = "$(FG)
	@echo " ar. option   = "$(AR_OPTION)
	@echo " MOD direc.   = "$(MOD_DIR)
	@echo " EXE direc.   = "$(EXE_DIR)
	@echo " LIB direc.   = "$(LIB_DIR)
	@echo " "

d_blas:
	(cd blas; make)
d_lpack:
	(cd lpack; make)
d_dust:
	(cd dust; make)
d_mpi:
	(cd mpi_output; make )
d_tools:
	(cd tools; make )
d_nm_subs:
	(cd new_main/subs; make )
d_plane:
	(cd plane; make )
d_stark:
	(cd stark; make )
d_subs:
	(cd subs; make )
d_unix:
	(cd unix; make )
d_newsubs:
	(cd newsubs; make ) 

# The following will create the executables

d_lte:
	(cd lte_hydro; make)
d_new_main:
	(cd new_main; make )
#
# To use the following command enter"
#         make -i clean
# The -i is necessary in case some files don't exist.

clean:
	rm -f lib/*.a
	rm -f mod/*.mod
	rm -f exe/*.exe
	rm -f */*.o
	rm -f */*/*.o
	rm -f */*/*/*.o
	rm -f */*.mod
	rm -f */*/*.mod
	rm -f */*/*/*.mod

dbg_clean:
	rm -f lib_dbg/*.a
	rm -f mod_dbg/*.mod
	rm -f exe_dbg//*exe
	rm -f */*.o
	rm -f */*/*.o
	rm -f */*/*/*.o
	rm -f */*.mod
	rm -f */*/*.mod
	rm -f */*/*/*.mod
