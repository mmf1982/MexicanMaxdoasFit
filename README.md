# MexicanMaxdoasFit
Inversion code for profile retrieval from slant column densities

You need vlidort2.7 (check here on how to receive it: http://www.rtslidort.com/mainprod_vlidort.html) to be able to run MexicanMaxdoasFit (MMF). 

## To compile the fortran core of MMF together with vlidort:
  * put the RETRIEVAL folder inside the vlidort folder with folders fo_main_1p4, util, vlidort_main, vlidort_def, vsup, vlidort_s_test
  * put the files makefile and make_code.py in that same folder.
  * create empty folders mod and obj inside that same folder
  * adjust some limits in the vlidort_pars.f90 file which is located inside the vlidort_def folder. Limits to change are e.g. MAXLAYERS (change from 26 to something around 100) and MAX_USER_VZANGLES (change to 35).
  * run python3.7 make_code.py from within that folder

## To run the retrieval of the included test measurements:
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
  * within the main folder, run: python3 -m pythoncode.MMFprofiling 

This produces a file netCDF4 test.nc which represents the retrieval of the data supplied here inside TESTDATA/MEAS.yml using the ancilary data (temperature, pressure, aerosol properties, etc.) given in TESTDATA/ANCIL.yml using the settings as specified in TESTDATA/SETTINGS.yml and inside the settings files that are linked inside that file ans also provided here in ancilMMF and configMMF. 

Please note that the interface to the example files provided here under TESTDATA is only an example. I chose this to be yml so that the data can be directly looked at as a plain text. I recommend a netCDF interface.

## Description of files
### configMMF/mmf_output.yml
This is a file describing the output variables inside the output file, something like a blueprint. Variables can be removed. Variables in the group FROM_ORIG need to be present in the original measurement dictionarry. The name inside the outputfile (given under the key "name") can be adapted", as well as the "description" or "long_name".

### configMMF/config.cnf 
This file describes output file settings, such as the output grid (not the retrieval grid), fillvalues and default values to use

### configMMF/mmf_settings.yml
flag limits and names for other settings. Retrieval settings are not set here, only settings for flags and furhter names of files

### ancilMMF/station/setup.txt
This is the blueprint setup file for the fortran core. The retrieval grid (and simulation grid) are set here. A proper description is yet to come. Some comments: Retrieval settings are provided here. For *numer of layers*, the first number is the total number of layers, the second is up to which number the retrieval is performed. *Sa read in?* currently only works for False. 

### ancilMMF/station/aerosol_apriori.txt, apriori_heights.txt, hcho_apriori.txt, no2_apriori.txt
 These are the apriori files to use. Note that the height grid is in meters and starts with the number of layers. The other files also start with the number of layers as the first entry. Further, the height grid is increasing. This grid is internally interpolated to the simulation grid specified in setup.txt. Note further that the aerosol a priori is internally rescaled to the provided aod value.
 
### TESTDATA
These contain examples of how ancillary information such as temperaure, pressure, aerosol properties, surface alebdo etc should be supplied (ANCIL.yml) as well as how the measurements should be supplied (MEAS.yml) and how the file structure etc. could be set up (SETTINGS.yml). The former two are only supplied here as yml files for ease of viewing. The corresponding wrapper for this is the supplied function "test_run" in MMFprofiling.py. The important part here is to pass three dicts and and output file name to mmf_stand_alone.
 
## suggested Versions
I suggest to use the following versions:

* python >= 3.7.5
* numpy >= 1.16.3
* yaml >= 5.1.2
* netCDF4 >= 1.5.1.2
* cf_units >= 2.1.3
* scipy >= 1.3.0
* gfortran >= 5.4.0

it might work with lower version as well.
