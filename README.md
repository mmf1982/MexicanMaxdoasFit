# MexicanMaxdoasFit
Inversion code for profile retrieval from slant column densities

You need vlidort2.7 to be able to run MexicanMaxdoasFit (MMF). 

To compile the fortran core of MMF together with vlidort:
  * put the RETRIEVAL folder inside the vlidort folder with folders fo_main_1p4, util, vlidort_main, vlidort_def, vsup, vlidort_s_test
  * put the files makefile and make_code.py in that same folder.
  * run python3.7 make_code.py from within that folder
    (you might need to adjust some limits in the vlidort_pars.f90 file, e.g. layer limit or angle limit)

To run the retrieval of the included test measurements:
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



   
