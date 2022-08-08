
#
# Define some variables
#

UTIL_PATH = util
VSUP_PATH = vsup
VLID_DEF_PATH = vlidort_def
VLID_MAIN_PATH = vlidort_main
VLID_TEST_PATH = vlidort_s_test
VLID_RET_PATH = RETRIEVAL
FO_MAIN_PATH = fo_main_1p4

MOD_PATH = mod
OBJ_PATH = obj

MOD_FILES = $(MOD_PATH)/*.mod
OBJ_FILES = $(OBJ_PATH)/*.o

#
# Define default shell make will use
#

SHELL = /bin/bash


#
# Define FORTRAN90 compiler to use (can be defined here, but usually defined on makefile command line)
#

#  Intel
# FC = ifort

#  gfortran
FC = gfortran

#  g95
#FC = g95

#  NAG
#FC = f95

#
# Define FORTRAN90 compiler flags
#

FFLAGS = -c

# Additional flags for Intel
ifeq ($(FC), ifort)
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -module $(MOD_PATH)
	FFLAGS_DEBUG = -g -warn all -debug -check all -traceback -O0 -check bounds -check uninit -ftrapuv -debug all 
	FFLAGS_OPENMP = -openmp
	FFLAGS_OPT = -O3
endif

# Additional flags for gfortran
ifeq ($(FC), gfortran)
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -J$(MOD_PATH) -mcmodel=large
	FFLAGS_DEBUG = -g -Wall -fbounds-check
	FFLAGS_OPT = -O3
#	FFLAGS_DEBUG = -g -Wall -fbounds-check -fbacktrace
endif

# Additional flags for g95
ifeq ($(FC), g95)
#      older g95
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -fmod=$(MOD_PATH)
	FFLAGS_DEBUG = -g -Wall -fbounds-check
#      g95 v0.92
#	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -fmod=$(MOD_PATH)
#	FFLAGS_DEBUG = -g -Wall -fbounds-check -ftrace=full
endif

# Additional flags for NAG
ifeq ($(FC), f95)
	FFLAGS := $(FFLAGS) -mdir $(MOD_PATH) -I$(MOD_PATH)
	FFLAGS_DEBUG = -w=obs -w=unused -C=array -C=undefined -gline
#	FFLAGS_DEBUG = -g -C=all -C=undefined -gline -mtrace=all -nan
#	FFLAGS_DEBUG = -w=obs -w=unused -g -C=all -C=undefined -gline -mtrace=all -nan
endif

# For debug build, use "make DEBUG=t"
ifeq ($(DEBUG), t)
	FFLAGS := $(FFLAGS) $(FFLAGS_DEBUG) -DDBG
endif

# For detailed output, i.e. convergence steps use "make DETAIL=t"
ifeq ($(DETAIL), t)
	FFLAGS := $(FFLAGS) -DDETAIL -DDBG
endif

# For optimized build, use "make OPT=t"
ifeq ($(OPT), t)
	FFLAGS := $(FFLAGS) $(FFLAGS_OPT)
endif

# To build TG retrieval code, use "TG=t"
ifeq ($(TG), t)
	FFLAGS := $(FFLAGS) -DTG
	MTEST="NO2_profile_wf.exe"
endif

# To build aerosol retrieval code, use "AEROSOL=t"
ifeq ($(AEROSOL),t)
	FFLAGS := $(FFLAGS) -DAEROSOL
	MTEST="AEROSOL_profile_wf.exe"
endif

# To build TG forward only, use "TGF=t"  !currently disabled
ifeq ($(TGF), t)
	FFLAGS := $(FFLAGS) -DTG
	FFLAGS := $(FFLAGS) -DFORW
	MTEST="makemokspectra.exe"
endif

# To build aerosol forward only, use "AEROSOLF=t"  !currently disabled
ifeq ($(AEROSOLF),t)
	FFLAGS := $(FFLAGS) -DAEROSOL
	FFLAGS := $(FFLAGS) -DFORW
	MTEST="makemokspectra_O4.exe"
endif

# To build retrieval code that uses retrieval in log space, use "LOG=t"
ifeq ($(LOG),t)
	FFLAGS := $(FFLAGS) -DLOGSPACE
endif

# To build TG that uses VMR as apriori, use "VMR=t".
ifeq ($(VMR),t)
	FFLAGS := $(FFLAGS) -DVMRIN
endif

# To build forward, use "FORW=t"
ifeq ($(FORW),t)
	FFLAGS := $(FFLAGS) -DFORW
	MTEST="forward.exe"
endif


.SUFFIXES:

#
# Define list of source files
# (Note: ordering is important because of dependencies)
#

BASE_SOURCES =
SOURCES =
L_SOURCES =
LPS_SOURCES =
LCS_SOURCES =
PROFILE_SOURCES =
N02_SOURCES =
AEROSOL_SOURCES =

BASE_SOURCES +=   \
   $(VLID_DEF_PATH)/vlidort_pars.f90

SOURCES +=   \
   $(BASE_SOURCES) \
   $(VLID_MAIN_PATH)/lapack_tools.f90		\
   $(VLID_DEF_PATH)/vlidort_inputs_def.f90	\
   $(VLID_DEF_PATH)/vlidort_sup_brdf_def.f90	\
   $(VLID_DEF_PATH)/vlidort_sup_ss_def.f90	\
   $(VLID_DEF_PATH)/vlidort_sup_sleave_def.f90	\
   $(VLID_DEF_PATH)/vlidort_sup_def.f90		\
   $(VLID_DEF_PATH)/vlidort_outputs_def.f90	\
   $(VLID_DEF_PATH)/vlidort_io_defs.f90		\
   $(VLID_DEF_PATH)/vlidort_work_def.f90	\
   $(VLID_MAIN_PATH)/vlidort_aux.f90		\
   $(VLID_MAIN_PATH)/vlidort_getplanck.f90	\
   $(VLID_MAIN_PATH)/vlidort_geometry.f90       \
   $(VLID_MAIN_PATH)/vlidort_Taylor.f90         \
   $(VLID_MAIN_PATH)/vlidort_inputs.f90		\
   $(VLID_MAIN_PATH)/vlidort_miscsetups.f90	\
   $(VLID_MAIN_PATH)/vlidort_multipliers.f90	\
   $(VLID_MAIN_PATH)/vlidort_corrections.f90	\
   $(VLID_MAIN_PATH)/vlidort_thermalsup.f90	\
   $(VLID_MAIN_PATH)/vlidort_solutions.f90	\
   $(VLID_MAIN_PATH)/vlidort_bvproblem.f90	\
   $(VLID_MAIN_PATH)/vlidort_intensity.f90	\
   $(VLID_MAIN_PATH)/vlidort_writemodules.f90	\
   $(VLID_MAIN_PATH)/vlidort_pack.f90		\
   $(VLID_MAIN_PATH)/vlidort_unpack.f90		\
   $(VLID_MAIN_PATH)/vlidort_masters.f90

L_SOURCES += \
   $(VLID_DEF_PATH)/vlidort_lin_inputs_def.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_sup_brdf_def.f90  \
   $(VLID_DEF_PATH)/vlidort_lin_sup_ss_def.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_sup_sleave_def.f90\
   $(VLID_DEF_PATH)/vlidort_lin_sup_def.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_outputs_def.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_io_defs.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_work_def.f90	  \
   $(VLID_MAIN_PATH)/vlidort_l_inputs.f90	  \
   $(VLID_MAIN_PATH)/vlidort_la_miscsetups.f90	  \
   $(VLID_MAIN_PATH)/vlidort_la_corrections.f90	  \
   $(VLID_MAIN_PATH)/vlidort_ls_corrections.f90   \
   $(VLID_MAIN_PATH)/vlidort_l_thermalsup.f90	  \
   $(VLID_MAIN_PATH)/vlidort_lpc_solutions.f90	  \
   $(VLID_MAIN_PATH)/vlidort_lpc_bvproblem.f90	  \
   $(VLID_MAIN_PATH)/vlidort_lbbf_jacobians.f90   \
   $(VLID_MAIN_PATH)/vlidort_ls_wfsurface.f90	  \
   $(VLID_MAIN_PATH)/vlidort_ls_wfsleave.f90	  \
   $(VLID_MAIN_PATH)/vlidort_l_writemodules.f90 \
   $(VLID_MAIN_PATH)/vlidort_l_pack.f90		  \
   $(VLID_MAIN_PATH)/vlidort_l_unpack.f90
   

LPS_SOURCES += \
   $(VLID_MAIN_PATH)/vlidort_lp_miscsetups.f90    \
   $(VLID_MAIN_PATH)/vlidort_lp_corrections.f90   \
   $(VLID_MAIN_PATH)/vlidort_lp_solutions.f90	  \
   $(VLID_MAIN_PATH)/vlidort_lp_bvproblem.f90     \
   $(VLID_MAIN_PATH)/vlidort_lp_wfatmos.f90       \
   $(VLID_MAIN_PATH)/vlidort_lp_pack.f90	  \
   $(VLID_MAIN_PATH)/vlidort_lp_unpack.f90	  \
   $(VLID_MAIN_PATH)/vlidort_lps_masters.f90

LCS_SOURCES += \
   $(VLID_MAIN_PATH)/vlidort_lc_miscsetups.f90    \
   $(VLID_MAIN_PATH)/vlidort_lc_corrections.f90   \
   $(VLID_MAIN_PATH)/vlidort_lc_solutions.f90	  \
   $(VLID_MAIN_PATH)/vlidort_lc_bvproblem.f90     \
   $(VLID_MAIN_PATH)/vlidort_lc_wfatmos.f90       \
   $(VLID_MAIN_PATH)/vlidort_lc_pack.f90	  \
   $(VLID_MAIN_PATH)/vlidort_lc_unpack.f90	  \
   $(VLID_MAIN_PATH)/vlidort_lcs_masters.f90

NO2_SOURCES += \
   $(VLID_RET_PATH)/constants.f90        \
   $(VLID_RET_PATH)/sorting.f90          \
   $(VLID_RET_PATH)/handle_setupfile.f90 \
   $(VLID_RET_PATH)/apriori_sa_read.f90 \
   $(VLID_RET_PATH)/stringmod.f90        \
   $(VLID_RET_PATH)/readinput.f90        \
   $(VLID_RET_PATH)/aerosolread.f90      \
   $(VLID_RET_PATH)/optpropin.F90        \
   $(VLID_RET_PATH)/atmosconditions.f90  \
   $(VLID_RET_PATH)/initialize.F90       \
   $(VLID_RET_PATH)/calcSCDK.f90         \
   $(VLID_RET_PATH)/gaussnewton.F90      \
   $(VLID_RET_PATH)/write_output.F90
   

AEROSOL_SOURCES += \
   $(VLID_RET_PATH)/constants.f90        \
   $(VLID_RET_PATH)/sorting.f90          \
   $(VLID_RET_PATH)/handle_setupfile.f90 \
   $(VLID_RET_PATH)/apriori_sa_read.f90 \
   $(VLID_RET_PATH)/stringmod.f90        \
   $(VLID_RET_PATH)/readinput.f90        \
   $(VLID_RET_PATH)/aerosolread.f90      \
   $(VLID_RET_PATH)/optpropin2.f90        \
   $(VLID_RET_PATH)/atmosconditions.f90  \
   $(VLID_RET_PATH)/initialize.f90       \
   $(VLID_RET_PATH)/calcSCDK.f90         \
   $(VLID_RET_PATH)/gaussnewton.F90      \
   $(VLID_RET_PATH)/write_output.F90


# (Include vector supplement source files)
include $(VSUP_PATH)/makefile.vsup


# this is copied from new makefile
include $(FO_MAIN_PATH)/makefile.fo


SOURCES_NO2_PROFILE = $(FO_SOURCES_L_Vector) + \
    $(SOURCES) +  \
   $(L_SOURCES) \
   $(LPS_SOURCES) \
   $(NO2_SOURCES) \
   $(VLID_RET_PATH)/writeinput.f90 \
   $(VLID_RET_PATH)/Profiling.F90

#NO2_profile_wf.f90
#Profile_WF.F90

SOURCES_FORWARD = $(FO_SOURCES_L_Vector) + \
    $(SOURCES) +  \
   $(L_SOURCES) \
   $(LPS_SOURCES) \
   $(NO2_SOURCES) \
   $(VLID_RET_PATH)/writeinput.f90 \
   $(VLID_RET_PATH)/Forward.F90

# these are currently deactivated
SOURCES_MAKEMOK =  $(FO_SOURCES_L_Vector) + \
    $(SOURCES) +  \
   $(L_SOURCES) \
   $(LPS_SOURCES) \
   $(NO2_SOURCES) \
   $(VLID_RET_PATH)/writeinput.f90 \
   $(VLID_RET_PATH)/makemokspectra.F90

# these are currently deactivated
SOURCES_TEST_O4 =  $(FO_SOURCES_L_Vector) + \
   $(SOURCES) +  \
   $(L_SOURCES) \
   $(LPS_SOURCES) \
   $(AEROSOL_SOURCES) \
   $(VLID_RET_PATH)/writeinput.f90 \
   $(VLID_RET_PATH)/O4_intensity_test.f90

# these are currently deactivated
SOURCES_MAKEMOK_O4 =  $(FO_SOURCES_L_Vector) + \
   $(SOURCES) +  \
   $(L_SOURCES) \
   $(LPS_SOURCES) \
   $(AEROSOL_SOURCES) \
   $(VLID_RET_PATH)/writeinput.f90 \
   $(VLID_RET_PATH)/makemokspectra_04.F90   


#
# Define pattern rules for creating object files:
#

#.SUFFIXES:

$(OBJ_PATH)/%.o : $(VLID_DEF_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VLID_MAIN_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VLID_TEST_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VLID_RET_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VLID_RET_PATH)/%.F90
	$(FC) $(FFLAGS) $< -o $@
# For utility source files
$(OBJ_PATH)/%.o : $(UTIL_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@


#
# Define object files
# notdir removes directory path
# filter removes white spaces
# patsubst 1, 2 replaces the ending 1 with ending 2
# addprefix 1, list  adds 1 infront of each element in list 

#F90SOURCES_SOLAR := $(notdir $(filter %.f90, $(SOURCES_SOLAR)))
#F90OBJECTS_SOLAR := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_SOLAR)))

#F90SOURCES_AEROSOL_PROFILE := $(notdir $(filter %.f90, $(SOURCES_AEROSOL_PROFILE)))
#F90OBJECTS_AEROSOL_PROFILE := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_AEROSOL_PROFILE)))

F90SOURCES_NO2_PROFILE := $(notdir $(filter %.f90 %.F90, $(SOURCES_NO2_PROFILE)))
F90SOURCES_NO2_PROFILE := $(patsubst %.F90, %.o,$(F90SOURCES_NO2_PROFILE))
F90OBJECTS_NO2_PROFILE := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_NO2_PROFILE)))

# currently not used
F90SOURCES_MAKEMOK     := $(notdir $(filter %.f90 %.F90, $(SOURCES_MAKEMOK)))
F90SOURCES_MAKEMOK     := $(patsubst %.F90, %.o,$(F90SOURCES_MAKEMOK))
F90OBJECTS_MAKEMOK     := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_MAKEMOK)))
# currently not used
F90SOURCES_MAKEMOK_O4  := $(notdir $(filter %.f90 %.F90, $(SOURCES_MAKEMOK_O4)))
F90SOURCES_MAKEMOK_O4  := $(patsubst %.F90, %.o, $(F90SOURCES_MAKEMOK_O4))
F90OBJECTS_MAKEMOK_O4  := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_MAKEMOK_O4)))
# currently not used
F90SOURCES_TEST_O4 := $(notdir $(filter %.f90, $(SOURCES_TEST_O4)))
F90OBJECTS_TEST_O4 := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_TEST_O4)))

F90SOURCES_FORWARD := $(notdir $(filter %.f90 %.F90, $(SOURCES_FORWARD)))
F90SOURCES_FORWARD := $(patsubst %.F90, %.o,$(F90SOURCES_FORWARD))
F90OBJECTS_FORWARD := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_FORWARD)))

#
# Define desired targets
#

#replace the two below by the next line
$(MTEST): $(F90OBJECTS_NO2_PROFILE)
	$(FC) $^ -o $@
#NO2_profile_wf.exe: $(F90OBJECTS_NO2_PROFILE)
#	$(FC) $^ -o $@
#AEROSOL_profile_wf.exe:  $(F90OBJECTS_AEROSOL_PROFILE)
#	$(FC) $^ -o $@

forward.exe: $(F90OBJECTS_FORWARD)
	$(FC) $^ -o $@

# currently not used
makemokspectra.exe: $(F90OBJECTS_MAKEMOK)
	$(FC) $^ -o $@
# currently not used
makemokspectra_O4.exe: $(F90OBJECTS_MAKEMOK_O4)
	$(FC) $^ -o $@
# currently not used
O4_intensity_test.exe: $(F90OBJECTS_TEST_O4)
	$(FC) $^ -o $@
	

.PHONY: clean
clean:
	rm -f *.o $(OBJ_FILES) *.mod $(MOD_FILES) *.log

