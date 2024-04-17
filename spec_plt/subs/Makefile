#Include file with system dependent definitions

include ../../Makefile_definitions

#
LIB=$(LIB_DIR)libspec.a
LIB2=$(LIB_DIR)libsubs.a
#
#
# We only make dispgen.exe if there has been a changen to the main *.f files,
# or if the libraries have changed since the last link.
#
##         $(LIB)(num_rec.o)\
#
OBJ    = $(LIB)(convolve.o)\
         $(LIB)(clip.o)\
         $(LIB)(compute_tau_radial.o)\
         $(LIB)(comp_ip_ew.o)\
         $(LIB)(create_gauss.o)\
         $(LIB)(divide_by_cont.o)\
         $(LIB)(do_int_rec_ap.o)\
         $(LIB)(do_sn_shell_shift.o)\
         $(LIB)(determ_reddening.o)\
         $(LIB)(fix_pop_oscillations.o)\
         $(LIB)(get_assoc_ion.o)\
         $(LIB)(get_mag.o)\
         $(LIB)(get_nd_nt_nit.o)\
         $(LIB)(hiabs.o)\
         $(LIB)(h2abs.o)\
         $(LIB)(int_rec_ap.o)\
         $(LIB)(is_line_abs_v2.o)\
         $(LIB)(rd_ion_locations.o)\
         $(LIB)(rd_ion_locations_v2.o)\
         $(LIB)(rd_sing_pol_i.o)\
         $(LIB)(read_corrections.o)\
         $(LIB)(set_ip_xspat_unit.o)\
         $(LIB)(uvabs_v2.o)\
         $(LIB)(adjust_jes_v2.o)\
         $(LIB)(filt_pass_band.o)\
         $(LIB)(get_prof_params_v2.o)\
         $(LIB)(cnvlv_es_one_par_v2.o)\
         $(LIB)(cnvlv_es_two_par_v2.o)\
         $(LIB)(cnvrt.o)\
         $(LIB)(cnvrt_j.o)\
         $(LIB)(dp_cnvrt.o)\
         $(LIB)(dp_cnvrt_j.o)\
         $(LIB)(dp_cnvrt_j_v2.o)\
         $(LIB)(plt_jh_opt_desc.o)\
         $(LIB)(raw_subphot_v2.o)\
         $(LIB)(rd_ew.o)\
         $(LIB)(rd_mod.o)\
         $(LIB)(recom_opac_v2.o)\
         $(LIB)(rd_row_data.o)\
         $(LIB)(rd_xy_data_usr.o)\
         $(LIB)(rd_monte_line.o)\
         $(LIB)(rite_2dj_asc_jh.o)\
         $(LIB)(dp_rd_obs_data_v2.o)\
         $(LIB)(rd_obs_data_v2.o)\
         $(LIB)(rd_zm_data.o)\
         $(LIB)(zm_data_module.o)
  
OBJ2    = $(LIB2)(poidev.o)\
          $(LIB2)(ran2.o)

all : $(LIB)
	cp $(INSTALL_DIR)perm_mod_file jnk.mod
	mv -f *.mod $(MOD_DIR)
	rm -f *.o

$(LIB) : $(OBJ) $(OBJ2)

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f
	$(F90) -c $(FG)  $<
	ar $(AR_OPTION) $(LIB) $*.o

$(LIB2)(%.o) : %.f
	$(F90) -c $(FG)  $<
	ar $(AR_OPTION) $(LIB2) $*.o

convolve.o : convolve.f parameters.inc constants.inc
num_rec.o  : num_rec.f   parameters.inc constants.inc
uvab.o     : uvabs.f     parameters.inc constants.inc

rd_xy_data.o : ../../subs/new_gen_in.o

clean :
	rm -f *.o

clean_lib :
	rm -f $(LIB)
