# ==============================================================================
# This file is part of GoMars since 2023.
#
# GoMars is a Martian general circulation model developed in Institute of
# Atmospheric Physics (IAP), Chinese Academy of Sciences (CAS).
#
# GMCORE is a dynamical core for atmospheric model used in GoMars.
#
# GoMars and GMCORE are distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY. You may contact authors for helping or cooperation.
# ==============================================================================

cmake_policy(SET CMP0076 NEW)

target_sources(gmcore PRIVATE
  gomars_v1_const_mod.F90
  gomars_v1_driver_mod.F90
  gomars_v1_namelist_mod.F90
  gomars_v1_objects_mod.F90
  gomars_v1_output_mod.F90
  gomars_v1_tracers_mod.F90
  gomars_v1_types_mod.F90
  gomars_v1_orbit_mod.F90
  gomars_v1_pbl_mod.F90
  gomars_v1_lsm_mod.F90
  gomars_v1_co2cyc_mod.F90
  gomars_v1_dustcyc_mod.F90
  gomars_v1_mp_mod.F90
  gomars_v1_cnvadj_mod.F90
  gomars_v1_damp_mod.F90
  gomars_v1_rad_mod.F90
  gomars_v1_utils_mod.F90
  nasa/sedim.F90
  nasa/microphys.F90
)