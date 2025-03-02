# Include local system definitions

include ../../Makefile_definitions

OPTION=$(FG)

# Library where object files will/are stored

LIB=$(LIB_DIR)libsubs.a
LIB_MOD=$(LIB_DIR)libmod.a

# Entries for library

#      $(LIB)(ba_ehb_ff_update_v1.o)\

OBJ = $(LIB)(add_par_to_full_mpi_v1.o)\
      $(LIB)(accel_tri_sol.o)\
      $(LIB)(adjust_cmat_to_lte.o)\
      $(LIB)(adjust_sn_r_grid.o)\
      $(LIB)(adjust_pops_mpi_v1.o)\
      $(LIB)(allocate_wse_arrays_mpi_v1.o)\
      $(LIB)(auto_clump_rev.o)\
      $(LIB)(auto_add_ion.o)\
      $(LIB)(ba_update_mpi_v1.o)\
      $(LIB)(check_ba_sol.o)\
      $(LIB)(check_eq_in_ba.o)\
      $(LIB)(check_hydro_def.o)\
      $(LIB)(check_ions_present.o)\
      $(LIB)(check_mpi_equality.o)\
      $(LIB)(check_tmin.o)\
      $(LIB)(chk_cool_rates_consist.o)\
      $(LIB)(cmf_blkband_mpi_v1.o)\
      $(LIB)(cmf_diag_band_mpi_v1.o)\
      $(LIB)(cmf_tri_band_mpi_v2.o)\
      $(LIB)(cmf_triband_thomas_mpi_v1.o)\
      $(LIB)(create_correction_sum_mpi_v1.o)\
      $(LIB)(create_iv_links_v2.o)\
      $(LIB)(def_new_rg_v1.o)\
      $(LIB)(determine_nse.o)\
      $(LIB)(det_lte_ed.o)\
      $(LIB)(ddt_work.o)\
      $(LIB)(do_cmf_hydro_mpi_v1.o)\
      $(LIB)(do_species_decays_v2.o)\
      $(LIB)(do_full_r_grid_v1.o)\
      $(LIB)(do_vel_regrid_v2.o)\
      $(LIB)(do_tau_regrid_v2.o)\
      $(LIB)(do_gam_abs_approx_v2.o)\
      $(LIB)(do_ng_accel_band_v4.o)\
      $(LIB)(do_wind_vel_v1.o)\
      $(LIB)(ba_ehb_bf_update_mpi_v1.o)\
      $(LIB)(ba_ehb_ff_update_mpi_v1.o)\
      $(LIB)(eval_adiabatic_mpi_v1.o)\
      $(LIB)(eval_rad_decay_mpi_v1.o)\
      $(LIB)(evalse_qwvj_mpi_v1.o)\
      $(LIB)(evalse_lowt_mpi_v1.o)\
      $(LIB)(evalse_x_qwvj_mpi_v1.o)\
      $(LIB)(eval_temp_ddt_mpi_v1.o)\
      $(LIB)(eval_shock_power.o)\
      $(LIB)(fiddle_pop_corrections_mpi_v1.o)\
      $(LIB)(fixpop_in_ba_v2.o)\
      $(LIB)(fixpop_t.o)\
      $(LIB)(fixpop_in_ba_v3.o)\
      $(LIB)(generate_full_matrix_v3.o)\
      $(LIB)(get_decay_lum_v1.o)\
      $(LIB)(get_edep_shock_power.o)\
      $(LIB)(get_excite_temp.o)\
      $(LIB)(get_non_loc_gamma_energy_v2.o)\
      $(LIB)(get_pops_at_prev_time_step_v5.o)\
      $(LIB)(old_grid_module.o)\
      $(LIB)(open_rw_eddfactor.o)\
      $(LIB)(out_line_force.o)\
      $(LIB)(prevent_low_t_mpi_v1.o)\
      $(LIB)(read_ba_data_mpi_v1.o)\
      $(LIB)(read_bcd_mat.o)\
      $(LIB)(read_time_model_v1.o)\
      $(LIB)(read_seq_time_file_v1.o)\
      $(LIB)(rv_sn_model_sniin.o)\
      $(LIB)(set_ba_storage_mpi_v1.o)\
      $(LIB)(set_ba_storage.o)\
      $(LIB)(set_depth_consistency.o)\
      $(LIB)(set_rv_hydro_model_v3.o)\
      $(LIB)(set_imp_vec.o)\
      $(LIB)(set_pops_to_pure_lte.o)\
      $(LIB)(set_sn_r_grid.o)\
      $(LIB)(set_new_grid_v2.o)\
      $(LIB)(smooth_pops_as_we_iterate.o)\
      $(LIB)(solveba_mpi_v1.o)\
      $(LIB)(spl_clump.o)\
      $(LIB)(steq_multi_mpi_v1.o)\
      $(LIB)(steqne_v4.o)\
      $(LIB)(steq_advec_mpi_v1.o)\
      $(LIB)(steq_co_mov_deriv_mpi_v1.o)\
      $(LIB)(store_old_grid.o)\
      $(LIB)(store_ba_data_mpi_v1.o)\
      $(LIB)(sub_guess_dc_v2.o)\
      $(LIB)(sum_steq_sol.o)\
      $(LIB)(tst_rd_seq.o)\
      $(LIB)(var_op_mpi_v1.o)\
      $(LIB)(var_x_opa_eta_mpi_v1.o)\
      $(LIB)(vehb_byj_v1.o)\
      $(LIB)(vsebyj_multi_mpi_v1.o)\
      $(LIB)(vsebyj_x_mpi_v1.o)\
      $(LIB)(write_cmat.o)\
      $(LIB)(write_seq_time_file_mpi_v1.o)\
      $(LIB)(update_ba_for_line_mpi_v1.o)\
      $(LIB)(wr_asci_steq_mpi_v1.o)\
      $(LIB)(wr_level_links.o)

#cmf_sub_sonic_hydro.o

all : $(LIB) XNUC XAUTO XCHG XTWO XNON
	cp -f $(INSTALL_DIR)perm_mod_file jnk.mod
	mv -f *.mod $(MOD_DIR)
	rm -f *.o

$(LIB) : $(OBJ)

STEQ=$(LIB_MOD)(steq_data_mod.o)
CHG=$(LIB_MOD)(chg_exch_mod_v3.o)
MOD=$(LIB_MOD)(mod_cmfgen.o)
DIS=$(LIB_MOD)(mod_lev_dis_blk.o)
CONT=$(LIB_MOD)(control_variable_mod.o)
NUC=$(LIB_MOD)(nuc_iso_mod.o)
HYD=$(LIB_MOD)(cmf_hydro_module.o)
EDD=$(LIB_MOD)(eddfac_rec_defs_mod.o)

UP_MODS= $(LIB_DIR)libmod.a(mod_cmfgen.o)\
       $(LIB_DIR)libmod.a(control_variable_mod.o)\
       $(LIB_DIR)libmod.a(line_mod.o)\
       $(LIB_DIR)libmod.a(line_vec_mod.o)\
       $(LIB_DIR)libmod.a(opac_mod.o)\
       $(LIB_DIR)libmod.a(radiation_mod.o)\
       $(LIB_DIR)libmod.a(steq_data_mod.o)\
       $(LIB_DIR)libmod.a(var_rad_mod.o)\
       $(LIB_DIR)libmod.a(wr_levels_links.o)


# Note: Not all files depend on mod_cmfgen and steq_data_mod, but assuming this dependency
# reduces signicantly the number of specific specifications.

$(LIB)(%.o) : %.f  $(STEQ) $(MOD) $(CONT)
	$(F90) -c $(OPTION) $<
	ar $(AR_OPTION) $(LIB) $*.o

# General rules to make obects in a library and object files
# Specific dependency of source files which uses external 
# fortran moudles. (The order of compilation matters.)

$(LIB)(adjust_pops_mpi_v1.o) : adjust_pops_mpi_v1.f $(LIB)(old_grid_module.o) $(MOD) $(CONT)
	$(F90) -c $(OPTION) adjust_pops_mpi_v1.f
	ar $(AR_OPTION) $(LIB) adjust_pops_mpi_v1.o

$(LIB)(determine_nse.o) : determine_nse.f $(STEQ) $(CHG) $(MOD)
	$(F90) -c $(OPTION) determine_nse.f
	ar $(AR_OPTION) $(LIB) determine_nse.o

$(LIB)(do_cmf_hydro_mpi_v1.o) : do_cmf_hydro_mpi_v1.f $(LIB)(old_grid_module.o) $(LIB)(cmf_hydro_module.o)
	$(F90) -c $(OPTION) do_cmf_hydro_mpi_v1.f
	ar $(AR_OPTION) $(LIB) do_cmf_hydro_mpi_v1.o

# Optimizaton seems not to work with PGF95 compiler for next routine

$(LIB)(eval_temp_ddt_v2.o) : eval_temp_ddt_v2.f $(STEQ) $(MOD)
	$(F90) -c $(OPTION)  eval_temp_ddt_v2.f
	ar $(AR_OPTION) $(LIB) eval_temp_ddt_v2.o

$(LIB)(open_rw_eddfactor.o) : open_rw_eddfactor.f $(EDD)
	$(F90) -c $(OPTION)  open_rw_eddfactor.f
	ar $(AR_OPTION) $(LIB) open_rw_eddfactor.o

$(LIB)(get_non_loc_gamma_energy.o) : get_non_loc_gamma_energy.f $(CONT) $(NUC)
	$(F90) -c $(OPTION) get_non_loc_gamma_energy.f
	ar $(AR_OPTION) $(LIB) get_non_loc_gamma_energy.o

$(LIB)(store_old_grid.o) : store_old_grid.f $(MOD)  $(LIB)(old_grid_module.o)
	$(F90) -c $(OPTION) store_old_grid.f
	ar $(AR_OPTION) $(LIB) store_old_grid.o

$(LIB)(update_ba_for_line.o) : update_ba_for_line.f $(UP_MODS)
	$(F90) -c $(OPTION)  update_ba_for_line.f
	ar $(AR_OPTION) $(LIB) update_ba_for_line.o

$(LIB)(var_op_v10.o) : var_op_mpi_v1.f $(DIS)
	$(F90) -c $(OPTION)  var_op_mpi_v1.f
	ar $(AR_OPTION) $(LIB) var_op_mpi_v1.o

$(LIB)(cmf_triband_thomas_mpi_v1.o) : cmf_triband_thomas_mpi_v1.f $(MOD)
	$(F90) -c $(OPTION)  cmf_triband_thomas_mpi_v1.f
	ar $(AR_OPTION) $(LIB) cmf_triband_thomas_mpi_v1.o

$(LIB)(wr_level_links.o) : wr_level_links.f $(MOD)
	$(F90) -c $(OPTION)  wr_level_links.f
	ar $(AR_OPTION) $(LIB) wr_level_links.o

$(STEQ):
	(cd ../mod_subs ; make)

$(MOD):
	(cd ../mod_subs ; make)

$(CONT):
	(cd ../mod_subs ; make)

$(EDD):
	(cd ../mod_subs ; make)

$(NUC):
	(cd ../mod_subs ; make)

$(UP_MODS):
	(cd ../mod_subs ; make)

$(DIS):
	(cd ../mod_subs ; make)

$(CHG):
	(cd chg ; make)

XNUC:
	(make -f Make_nuc)

XAUTO:
	(cd auto; make)
XCHG:
	(cd chg; make)
XTWO:
	(cd two; make)
XNON:
	(cd non_therm; make)

clean :
	rm -f *.o 

clean_lib :
	rm $(LIB)
