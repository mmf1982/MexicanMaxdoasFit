# MexicanMaxdoasFit
Inversion code for profile retrieval from slant column densities

You need vlidort2.7 to be able to run MexicanMaxdoasFit (MMF). 

# To compile the fortran core of MMF together with vlidort:
  * put the RETRIEVAL folder inside the vlidort folder with folders fo_main_1p4, util, vlidort_main, vlidort_def, vsup, vlidort_s_test
  * put the files makefile and make_code.py in that same folder.
  * run python3.7 make_code.py from within that folder
    (you might need to adjust some limits in the vlidort_pars.f90 file, e.g. layer limit or angle limit)

# To run the retrieval of the included test measurements:
  * create the follwing folder structure (i.e. copy the produced executables inside the main path):
      ```
      ancilMMF:
          2p6_VLIDORT_ReadInput.cfg
          HCHO_223.txt
          O4_298.txt
          no2_298K_vanDaele.xs
          station:
              aerosol_apriori.txt
              apriori_heights.txt
              hcho_apriori.txt
              no2_apriori.txt
              setup.txt
      configMMF:
          config.cnf
          mmf_output.yml
          mmf_settings.yml
      pythoncode:
          MMFprofiling.py
          Tools_MMF.py
          yaml_netcdf4.py
      TESTDATA:
          ANCIL.yml
          MEAS.yml
          SETTINGS.yml
      AEROSOL_profile_wf.exe
      NO2_profile_wf.exe
      ```
  * within the main folder, run: python3.7 -m pythoncode.MMFprofiling

This produces a file netCDF4 test.nc which represents the retrieval of the data supplied here inside TESTDATA/MEAS.yml using the ancilary data (temperature, pressure, aerosol properties, etc.) given in TESTDATA/ANCIL.yml using the settings as specified in TESTDATA/SETTINGS.yml and inside the settings files that are linked inside that file ans also provided here in ancilMMF and configMMF. 

# Description of files
## mmf_output.yml
This is a file describing the output variables inside the output file. Variables can be removed. Variables in the group FROM_ORIG need to be present in the original measurement dictionarry. The name inside the outputfile (given under the key "name") can be adapted", as well as the "description"

## config.cnf 
This file describes output file settings, such as the output grid (not the retrieval grid), fillvalues and default values to use

## mmf_settings.yml
flag limits and names for other settings. Retrieval settings are not set here, only settings for flags and furhter names of files

## setup.txt
This is the blueprint setup file for the fortran core. The retrieval grid (and simulation grid) are set here. A proper description is yet to come. Retrieval settings are provided here 

## aerosol_apriori.txt, apriori_heights.txt, hcho_apriori.txt, no2_apriori.txt
 These are the apriori files to use. Note that the height grid is in meters and starts with the number of layers. The other files also start with the number of layers as the first entry. Further, the heihgt grid is increasing. This grid is internally interpolated on the simulation grid specified in setup.txt. Note further that the aerosol a priori is internally rescaled to the provided aod value.
 

   
