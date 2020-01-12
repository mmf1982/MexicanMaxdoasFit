'''
run python3 make_code.py to compile both aerosol and tracegas executables
'''
import os

filelist = ["mod/optpropin.mod",
            "mod/atmosconditions.mod",
            "mod/calcscdk.mod ",
            "mod/constants.mod",
            "mod/gaussnewton.mod",
            "mod/handle_setupfile.mod",
            "mod/initialize.mod",
            "mod/no2read.mod",
            "mod/readinput.mod",
            "mod/sorting.mod",
            "mod/strings.mod",
            "mod/aerosolread.mod",
            "mod/precision.mod",
            "mod/apriori_sa_read.mod",
            "mod/write_output.mod",
            "mod/writeinput.mod",
            "obj/optpropin.o",
            "obj/atmosconditions.o",
            "obj/calcSCDK.o ",
            "obj/constants.o",
            "obj/gaussnewton.o",
            "obj/handle_setupfile.o",
            "obj/initialize.o",
            "obj/NO2read.o",
            "obj/readinput.o",
            "obj/sorting.o",
            "obj/stringmod.o",
            "obj/aerosolread.o",
            "obj/precmod.o",
            "obj/Profile_WF.o" ,
            "obj/apriori_sa_read.o",
            "obj/write_output.o",
            "obj/writeinput.o",
            "obj/Profiling.o",
            "AEROSOL_profile_wf.exe",
            "NO2_profile_wf.exe"
            ]

for mfile in filelist:
    command = "rm "+ mfile
    os.system(command)

os.system("make OPT=t TG=t LOG=t") # I want this
#os.system("make OPT=t TG=t DETAIL=t")


# os.system("make DEBUG=t TG=t LOG=t DETAIL=t")

for mfile in filelist[:-1]:
    command = "rm "+ mfile
    os.system(command)


os.system("make OPT=t AEROSOL=t LOG=t")  # I want this

#os.system("make OPT=t AEROSOL=t DETAIL=t")

#os.system("make DEBUG=t AEROSOL=t LOG=t DETAIL=t")

#for mfile in filelist[:-1]:
#    command = "rm "+ mfile
#    os.system(command)#
# to make the forward:
# os.system("make FORW=t ")

#os.system("make OPT=t AEROSOLF=t ")



