'''
Module to run fortran programs for MMF retrieval.

These were designed originally for work within FRM4DOAS and build upon the
routines developed for UNAM between 2014 and 2017.

Author:  Martina M Friedrich
Date:  August 2018
Version:  2019.09

MMF is a retrieval code for MAX-DOAS measurements. It uses VLIDORT as forward
model and uses a 2 step approach: First, o4 SCD are used for aerosol profile
retrieval. Next, this aerosol retrieval is fixed in the following trace gas
retrieval.
Both retrievals use a optimal estimation approach where the regulatization
matrix is constructed currently within the code from the a-priori profile.
Alternatively, a regulatization matrix can be supplied. However, this is not
yet tested for retrieval in logarithmic space.
Both steps use a Levenberg Marquad iteration scheme. On compile time, it is
decided wether the gas profile input is in VMR (in PPM) or in average
concentration at level center in molec/ cm**3.

The core MMF code takes plain text files as input and outputs plain text files.
It works on a scan-by-scan basis. Each time the algorithm is called, it needs
the following files as input:
setup.txt
temppress_data.txt
apriori_heights.txt
gas_apriori.txt aerosol_apriori.txt
qdoas.asc output file
'''

import os
import sys
import time
import copy
import shutil
from collections import OrderedDict
from multiprocessing import Process
import multiprocessing
from scipy.interpolate import interp2d
import yaml
import numpy as np

try:
    import Tools_MMF as TM
except ImportError:
    from . import Tools_MMF as TM

try:
    import yaml_netcdf4 as yn
except ImportError:
    from . import yaml_netcdf4 as yn

SETDICT = {}

def initialize(mfile):
    '''
    Initialize Tools_MMF with parameters in provided file, e.g. if debug

    Parameters:
    -----------
    mfile: string
        path to configuration file containing general configurations
        frm4doas_tropo.cnf:
        need at least:
            GSETTINGS["DEBUG"]
            REALDIMS["MAXAXIS"]
            SYSERRFRAC[gasname]
            DEFAULTS["aod_apriori"]
    '''
    _ = TM.initialize(mfile)
    return


def processwrapper(inputs, outputs):
    '''
    wrapper of calling any function for running on multiple processors.

    In this setup, it is not possible to give any arguments to the function
    call. Therefore, it might be necessary to have the function to be a
    class function and have access to instance specific values.

    Parameters:
    ----------
    inputs: multiprocessing.manager().queue object
        Each entry in queue is a tuple of (job to do, label or scannumber)
    outputs: multiprocessing.manager().queue object
        in the case of MMF, each entry in the queue is a tuple of
        (scan number, list)
        list is a list of:
            scd, avk, profile (for tracegas) or
            aod, scd, avk, profile (for aerosol)
    '''
    for func, pos in iter(inputs.get, "STOP"):
        if TM.GSETTINGS["DEBUG"]:
            TM.dwrite("proc wrapper "+str(pos))
        result = func()
        outputs.put(result)
    return


def clean_mmf(mdirectory):
    '''
    clean the mmf working directory of all files and sub directories

    Parameters:
    ----------
    mdirectory: string
        path to directory to remove
    '''
    shutil.rmtree(mdirectory)
    return


class MMFSetupFile(object):
    '''
    class for MMF setupfile handling

    The caracter indicating a header in that file is "!". The caracter for a
    new input key is "*".
    This file is processed by the fortran MMF implementation and hence the
    format cannot be changed easily. The "keys" have specific names, too.
    '''
    def __init__(self, filename):
        '''
        constructor of MMFSetupFile

        Parameters:
        ----------
        filename: string
            setup file to open incuding path

        class variables:
        ----------------
        content: list
            a line by line list of the file content
        mdict: dictionarry
            lines starting with * are interpreted as keys, following lines
            as the entry
        header: list of strings string
            original header that will be re-produced in the new file
        '''
        mfile = open(filename, "r")
        self.content = mfile.readlines()
        mfile.close()
        self.mdict = OrderedDict()
        self.process_file()
        self.header = [line for line in self.content if line[0] in "!"]

    def process_file(self):
        '''
        read input file into dictionarry.

        Lines that start with * are treated as the key, following lines, until
        the next line with * are treated as the value. This format is needed
        by the MMF fortran implementation.
        '''
        for line in self.content:
            if "*" in line[:3]:
                mkey = line
                self.mdict[mkey] = []
            else:
                try:
                    self.mdict[mkey].append(line)
                except UnboundLocalError:
                    pass
        for mkey in self.mdict.keys():
            self.mdict[mkey] = "".join(self.mdict[mkey])
        return

    def update(self, look_up_name, mentry):
        '''
        add entry to inputfile dictionarry, or update existing.

        Parameters:
        ----------
        look_up_name: string
            key for the dictionarry. Includes * and new line caracter \n
        mentry: float, string, list or array
            entry that should replace the entry under key string. If the
            key was not existent before, it is added.
        '''
        if isinstance(mentry, (list, np.ndarray)):
            mentry = "\n".join([str(ii) for ii in mentry])
        self.mdict[look_up_name] = str(mentry) + "\n"
        return

    def get_entry(self, look_up_name):
        '''
        extract entry from dictionarry

        Parameters:
        ----------
        look_up_name: string
            key for the dictionarry. Includes * and new line caracter \n

        Returns:
        --------
        mentry: string
            entry that should replace the entry under key string. If the
            key was not existent before, this crashes with KeyError
        '''
        mentry = self.mdict[look_up_name].lstrip("\n")
        try:
            mentry = float(mentry)
        except ValueError:
            pass
        return mentry

    def writefile(self, newfilename):
        '''
        write the inputfile dictionarry in the form of an input file

        Parameter:
        ----------
        newfilename: string
            filename including path to file to save
        '''
        mfile = open(newfilename, "w")
        for line in self.header:
            mfile.write(line)
        for mkey in self.mdict.keys():
            mfile.write(mkey)
            if isinstance(self.mdict[mkey], (list, np.ndarray)):
                _ = [mfile.write(line) for line in self.mdict[mkey]]
            else:
                mfile.write(str(self.mdict[mkey]))
        mfile.close()
        return


class MMFresults(dict):
    '''
    class to hold MMF results from either tg or aer, derived from dict

    MMFresults is a dictionarry with some additional functionality, such as
    to gather several entries (from different scans) in the values part into an
    array with a single key: from a list of dictionarries, it makes a
    dictionarry containing arrays.
    '''
    allkeys = ["aod", "aod_err", "vcd_err", "p_aod", "scd_sim", "scd_meas",
               "scd_meas_err", "profile", "apriori", "avk",
               "noiseerror", "smootherror", "combierror", "rms_dscd",
               "rms_r_dscd", "scan_number", "layerthickness", "air_column",
               "pressure", "temperature", "middleheights", "version",
               "ssa", "asy", "avk_col", "sm_err_cov", "n_err_cov"]
    # oneD are only those 1D keys that are in ae and tg retrieval
    oneD = ["dof", "vcd", "flag", "counter", "Sa_scaling", "Sa_corrlen"]

    def __init__(self):
        '''
        Constructor of MMFresults

        In initializes the dictionarry-type MMFresults with empty lists for
        a range of keys, fitting both aer and tg retrievals.
        '''
        super(MMFresults, self).__init__()
        for key in self.allkeys+self.oneD:
            self[key] = []

    def join(self, mlist):
        '''
        makes an array for each key of each entry in mlist.

        mlist is a list of MMFresults objects. This function extracts each key
        from this list and appends those to the value of the corresponding key
        in the calling object.

        Parameters:
        -----------
        mlist: list
            list of MMFresults, one for each scan. Ordered

        Returns:
        --------
        local_erc: integer
            a local error code set to MMFTG and MMFAE failed, if no "version".
        '''
        local_erc = 0
        for entry in mlist:
            for key in entry.keys():
                if np.array(entry[key]).size:  # false for []
                    self[key].append(entry[key])
        self["scan_number"] = np.array(self["scan_number"])
        morder = np.argsort(self["scan_number"])
        poplist = []
        for key in self.keys():
            if len(self[key]) > 0:
                self[key] = np.array(self[key])[morder]
            else:
                poplist.append(key)
        for key in poplist:
            self.pop(key)
        for key in self.keys():
            if key in ["layerthickness", "middleheights"]:
                # if key in "middleheights":
                # internally, units are km
                self[key] = np.nanmean(self[key], axis=0)*TM.CONV["KM2M"]
        if "version" in self.keys():
            try:
                self["version"] = next(
                    (entr for entr in self["version"] if entr is not None))
            except:  # this means none of the ae of mmf have passed.
                self["version"] = np.nan
                local_erc = 1
        return local_erc


class MMF(object):
    '''
    Main class to run the MMF fortran routines for a single scan

    It writes the input asc file, the setup file, the tp file and can run
    tracegas retrieval or aerosol retrieval
    '''
    def __init__(self, workdir, tph, instr_height, scannumber, mmf_spec,
                 gasdict, raa, sza, el_ang, ae_profile, ae_tau, xgas,
                 mmfinputs, tgexe, aeexe, ae_asym, ae_ssa, albedo):
        '''
        Parameters:
        -----------
        workdir: string
            path to the working directory where the fortran program writes to,
            a unique part will be added
        tph: 2D array
            each lines has entries: altitude (km), pressure (hPa),  temp (K)
        instr_height: string
            instrument height in km
        scannumber: integer
            scan number
        mmf_spec: string
            path to mmf station specific input, such as: setup blueprint and
            a-priori files.
        gasdict: dictionarry
            having at least members: dscd, dscd_err, wavelength, name
        raa: float
            relative solar zenith angle (average for scan)
        sza: float
            solar zenith angle (average for scan)
        el_ang: 1D array
            elevation angle sequence, can be masked array
        ae_profile: 1D array or None
            if ae retrieval, None, if tg retrieval the ae profile on grid in
            setup file this needs to be the same as the ae retrieval
        ae_tau: float
            if ae retrieval, this is default value, else retrieved
        xgas: str
            internal name of gas, <gas>_<wavelength>nm
        mmfinputs: str
            path to directory with general MMF inputs (vlidort cfg & xs files)
        tgexe: str
            path to tg retrieval executable
        aeexe: str
            path to ae retrieval executable
        ae_asym: scalar float
            aerosol asymmetry parameter
        ae_ssa: scalar float
            aerosol single scattering albedo
        albedo: scalar float
            surface albedo
        '''
        try:
            os.makedirs(workdir, mode=0o775)
        except FileExistsError:
            pass
        self.workdir = workdir + str(scannumber) + "_" + xgas + "/"
        try:
            os.makedirs(self.workdir, mode=0o775)
        except FileExistsError:
            pass
        self.numang = TM.REALDIMS["MAXAXIS"]
        self.scan = scannumber
        self.errm = self.writeasc(gasdict, raa, sza, el_ang)
        self.writetph(tph)
        self.ae_tau = ae_tau
        self.mmfinputs = mmfinputs
        self.tgex = tgexe
        self.aeex = aeexe
        self.mmfinputs_specific = mmf_spec
        self.ssa = ae_ssa
        self.asy = ae_asym
        interm = self.writesetup(gasdict["wavelength"], instr_height,
                                 ae_profile, ae_tau, xgas, albedo)
        self.confname = interm[0]
        self.sa_scaling = interm[1]
        self.sa_correlation_length = interm[2]
        self.totgrid_layers = interm[3]
        self.systerrfraction = TM.SYSERRFRAC[gasdict["name"]]

    def writetph(self, tph):
        '''
        write temperature pressure height as txt file for fortran code

        Parameters:
        -----------
        tph: 2D array
            one line per altitude level, each line: altitude/ km, P/ hPa, T/ K
        '''
        with open(self.workdir + SETDICT["settings"]["TPNAME"], 'w') as fid:
            fid.write("altitude (km) \t pressure (hPa) \t temp (K) \n")
            for line in tph:
                fid.write("  ".join(line) + "\n")
        return

    def writesetup(self, wavl, instr_height, ae_profile, ae_tau, xgas, albedo):
        '''
        write setup txt file for fortran code

        Parameters:
        -----------
        wavl: 1x1 array
            wavelength in nm
        instr_height: str
            str of instrument height in km
        ae_profile: None or 1D array of ae_profile on grid as in setup (last x
            layers as given as 2nd number under * number of layers)
        ae_tau: float
            aerosol total tau. if ae retrieval, this is a priori aod
        xgas: str
            name of gas <gas>_<wavelength>nm
        albedo: 1x1 array
            surface albedo

        Returns:
        --------
        confname: str
            path to MMF fortran setup file
        sa_scale: float
            factor of a priori on the Sa
        sa_corr_len: float
            correlation length for Sa construction in km
        rlayers: int
            number of retrieval layers
        '''
        ascname = self.workdir + "/QDOAS.asc"
        confname = self.mmfinputs_specific + SETDICT["settings"]["MMFSETUP"]
        msetup = MMFSetupFile(confname)
        msetup.update("* input filename\n", '"'+ascname+'"')
        msetup.update("* surface albedo\n", str(np.squeeze(albedo)))
        msetup.update("* output directory\n", '"' + self.workdir + '"')
        msetup.update("* temperature pressure filename\n",
                      '"' + self.workdir + SETDICT["settings"]["TPNAME"] + '"')
        msetup.update("* wavelength in nm\n", wavl)
        msetup.update("* altitude in km\n", instr_height)
        msetup.update("* cross section folder\n",
                      '"' + self.mmfinputs + SETDICT["SIG_DIC"][
                          xgas.split("_")[0].lower()]+'"')
        sa_scale = msetup.get_entry("* scaling parameter for Sa matrix\n")
        sa_corr_len = msetup.get_entry("* Sa correlation length in km\n")
        msetup.update("* heights for apriori profile filename\n",
                      '"' + self.mmfinputs_specific + 'apriori_heights.txt"')
        msetup.update("* vlidort config filename\n",
                      '"' + self.mmfinputs + '2p6_VLIDORT_ReadInput.cfg"')
        rlayers = int(msetup.get_entry("* number of layers\n").split("\n")[1])
        if ae_profile is not None:
            msetup.update("* aerosol profile shape\n", ae_profile)
            msetup.update("* aerosol ext t, unisotropy g, omega w\n",
                          np.array([str(np.squeeze(ae_tau)),
                                    str(np.squeeze(self.asy)),
                                    str(np.squeeze(self.ssa))]))
            if "hcho" in xgas.lower():
                msetup.update("* apriori profile filename\n",
                              '"'+self.mmfinputs_specific + 'hcho_apriori.txt"')
            elif "no2" in xgas.lower():
                msetup.update("* apriori profile filename\n",
                              '"'+self.mmfinputs_specific + 'no2_apriori.txt"')
            else:
                msetup.update("* apriori profile filename\n",
                              '"'+self.mmfinputs_specific + 'gas_apriori.txt"')
        else:
            msetup.update("* aerosol ext t, unisotropy g, omega w\n",
                          np.array(
                              [str(np.squeeze(ae_tau)),
                               str(np.squeeze(self.asy)),
                               str(np.squeeze(self.ssa))]))
            msetup.update("* apriori profile filename\n",
                          '"' + self.mmfinputs_specific +
                          'aerosol_apriori.txt"')
        confname = self.workdir + "/" + SETDICT["settings"]["MMFSETUP"]
        msetup.writefile(confname)
        return confname, sa_scale, sa_corr_len, rlayers

    def writeasc(self, gasdict, raa, sza, el_ang):
        '''
        write plain text dscd file for fortran code

        Parameters:
        -----------
        gasdict: dict
            dict with keys dscd, dscd_err, wavlength and name
        raa: float
            relative azimuth angle average for scan in degree
        sza: float
            solar zenith angle average for scan in degree
        el_ang: 1D array
            elevation anlge sequence, starting with 2 zenith measurements

        Returns:
        --------
        mreturn: string if error, otherwise None

        Comments:
        ---------
        The format needed for the fortran routine in the QDOAS.asc file is:
        integer of number of lines (number of el ang measurements) followed by:
            sza, raa, el ang, dscd, dscd err
        a minimum dscd err of 0.5 percent of the dscd is assumed
        '''
        ascname = self.workdir + "QDOAS.asc"
        # mask = ~np.isnan(gasobj.dscd)
        # mask = np.isnan(gasobj.dscd)
        # I am currently not using any mask here
        if sza >= 90.0 or sza <= 0.0:
            mreturn = "[MMFprofiling.py] ERROR, sza is out of range: "+str(sza)
        else:
            mreturn = None
        raa_matrix = np.ones(el_ang.shape)*raa
        sza_matrix = np.ones(el_ang.shape)*sza
        if isinstance(gasdict["dscd"], np.ma.MaskedArray):
            dscd = gasdict["dscd"].filled(np.nan)
        else:
            dscd = gasdict["dscd"]
        if isinstance(gasdict["dscd_err"], np.ma.MaskedArray):
            dscd_err = gasdict["dscd_err"].filled(np.nan)
        else:
            dscd_err = gasdict["dscd_err"]
        # assume the error is never ever smaller than 0.1% of the measured dscd
        dscd_err = np.fmax(dscd_err, 0.001*dscd)
        matrix = np.stack(
            (sza_matrix, raa_matrix, el_ang, dscd, dscd_err), axis=1)
        notnan = matrix.shape[0]
        # -2  # now exclude just the two first angles which are zenith
        np.savetxt(ascname, matrix, header=str(notnan), delimiter=",",
                   comments="")
        return mreturn

    def run_ae(self):
        '''
        run the aerosol retrieval fortran code and load the plain text files

        Returns:
        --------
        results: dict
            dict containing ae results for a single scan, keys (see MMFresults)
            aod, aod_err, vcd_err, p_aod, scd_sim, scd_meas_err, profile,
            apriori, avk, noiseerror, smootherror, combierror, rms_dscd,
            rms_r_dscd, scan_number, layerthickness, air_column, pressure,
            temperature, middleheights, version, ssa, asy, avk_col, sm_err_cov,
            n_err_cov, dof, vcd, flag, counter, Sa_scaling, Sa_corrlen

        TODO:
            this and run_tg should be put together and called with a key, just
            differences should be condisdered
        '''
        results = MMFresults()
        results["version"] = None
        if self.errm is not None:
            TM.dwrite("aerosol " + str(self.errm)+" in scan " + str(self.scan))
            results["aod"] = np.nan
            results["aod_err"] = np.nan
            results["p_aod"] = np.nan
            scd = np.ones([self.numang, 4])*np.nan
            results["avk"] = np.ones(
                [self.totgrid_layers, self.totgrid_layers])*np.nan
            results["sm_err_cov"] = np.copy(results["avk"])
            results["n_err_cov"] = np.copy(results["avk"])
            profile = np.ones([self.totgrid_layers, 7])*np.nan
            results["noiseerror"] = np.ones([self.totgrid_layers])*np.nan
            results["smootherror"] = np.ones([self.totgrid_layers])*np.nan
            results["combierror"] = np.ones([self.totgrid_layers])*np.nan
            for key in results.oneD:
                results[key] = np.nan
            results["avk_col"] = np.nan * np.ones(self.totgrid_layers)
        else:
            try:
                os.system(
                    self.aeex + " " + self.workdir + " " + self.confname +
                    " " + str(self.scan))
                for key in ["version", "aod", "vcd", "dof", "flag", "counter"]:
                    results[key] = np.loadtxt(self.workdir+"/"+key+".dat")
                profile = np.loadtxt(self.workdir + '/height_profiles.dat',
                                     skiprows=2)[:, [0, 1, 4, 5, 6, 7, 8]]
                # I keep: pAOD, layerheight, air column, P, T, height a.s.l
                for key in ["noiseerror", "smootherror", "avk", "sm_err_cov",
                            "n_err_cov"]:
                    results[key] = np.loadtxt(
                        self.workdir+"/"+key+".dat", skiprows=2)[:]
                try:
                    scd = np.loadtxt(
                        self.workdir + '/SCD_table.dat', skiprows=2)[:, :4]
                except IndexError:  # only one off-axis angle was used
                    scd = np.loadtxt(
                        self.workdir + '/SCD_table.dat', skiprows=2)[:4]
                    scd = scd.reshape([1, scd.size])
                    TM.write_error("only one off-axis" +
                                   str(scd.shape), exc="", error=False)
                results["Sa_scaling"] = self.sa_scaling
                results["Sa_corrlen"] = self.sa_correlation_length
                results["avk_col"] = results["avk"].sum(axis=1)
                results["aod_err"] = self.make_total_err(
                    results["sm_err_cov"], results["n_err_cov"],
                    results["aod"])
            except OSError as erro:  # this means, mmf failed for some reason.
                if TM.GSETTINGS["DEBUG"]:
                    TM.write_error("MMF failed ae in " + str(self.workdir),
                                   erro, False, False)
                results["aod"] = np.nan
                results["aod_err"] = np.nan
                results["p_aod"] = np.nan
                scd = np.ones([self.numang, 4])*np.nan
                results["avk"] = np.ones(
                    [self.totgrid_layers, self.totgrid_layers])*np.nan
                results["sm_err_cov"] = np.copy(results["avk"])
                results["n_err_cov"] = np.copy(results["avk"])
                profile = np.ones([self.totgrid_layers, 7])*np.nan
                results["noiseerror"] = np.ones([self.totgrid_layers])*np.nan
                results["smootherror"] = np.ones([self.totgrid_layers])*np.nan
                results["combierror"] = np.ones([self.totgrid_layers])*np.nan
                for key in results.oneD:
                    results[key] = np.nan
                results["avk_col"] = np.nan * np.ones(self.totgrid_layers)
        results["apriori"] = profile[:, 1]
        results["layerthickness"] = profile[:, 2]
        results["air_column"] = profile[:, 3]
        results["pressure"] = profile[:, 4]*TM.CONV["Pa2hPa"]  # internal is Pa
        results["temperature"] = profile[:, 5]
        results["middleheights"] = profile[:, 6]
        results["profile"] = profile[:, 0]
        results["combierror"] = np.sqrt(
            results["noiseerror"]**2+results["smootherror"]**2 +
            (results["profile"]*self.systerrfraction)**2)
        if not np.isnan(results["vcd"]):
            lowest = np.argmin(results["middleheights"])
            lowest = results["middleheights"][lowest]-0.5*results[
                "layerthickness"][lowest]
            which = results["middleheights"] - lowest < 4.0
            results["p_aod"] = np.nansum(
                profile[which, 0] *
                profile[which, 2])
        results["el_ang"] = scd[:, 0]
        results["scd_sim"] = scd[:, 1]
        results["scd_meas"] = scd[:, 2]
        results["scd_meas_err"] = scd[:, 3]
        results["ssa"] = self.ssa
        results["asy"] = self.asy
        results["scan_number"] = self.scan
        if np.isnan(results["aod"]):
            results["flag"] = SETDICT["flags"]["fail"]
            results["rms_dscd"] = np.nan
            results["rms_r_dscd"] = np.nan
        else:
            results["flag"], results["rms_dscd"], results["rms_r_dscd"] = self.new_make_flag(
                results["flag"], results["scd_sim"], results["scd_meas"],
                results["scd_meas_err"], results["dof"], results["el_ang"])
        if TM.GSETTINGS["DEBUG"]:
            TM.dwrite("doing scan " + str(self.scan))
        clean_mmf(self.workdir)
        del results["el_ang"]
        return results

    def run_tg(self):
        '''
        run the trace gas retrieval fortran code and load the plain text files

        Returns:
        --------
        results: dict
            results dictionarry with the same keys as run_ae, but not all have
            content.

        TODO:
            see run_ae TODO
        '''
        results = MMFresults()
        if self.errm is not None:
            TM.write_error("MMF failed tg in scan "+str(self.scan), "", False)
            results["apriori"] = np.ones(self.totgrid_layers)*np.nan
            results["profile"] = np.ones(self.totgrid_layers)*np.nan
            results["noiseerror"] = results["profile"]
            results["smootherror"] = results["profile"]
            results["combierror"] = results["profile"]
            scd = np.ones([self.numang, 4])*np.nan
            results["avk"] = np.ones(
                [self.totgrid_layers, self.totgrid_layers])*np.nan
            results["sm_err_cov"] = np.copy(results["avk"])
            results["n_err_cov"] = np.copy(results["avk"])
            for key in results.oneD:
                results[key] = np.nan
            results["vcd_err"] = np.nan
            results["avk_col"] = np.nan * np.ones(self.totgrid_layers)
        else:
            try:
                if np.isnan(self.ae_tau):
                    raise ValueError("tau is nan")
                os.system(self.tgex + " " + self.workdir + " "+self.confname +
                          " " + str(self.scan))
                profile = np.loadtxt(self.workdir + '/height_profiles.dat',
                                     skiprows=2, usecols=[1, 4])
                results["apriori"] = profile[:, 0]
                results["profile"] = profile[:, 1]
                for key in ["noiseerror", "smootherror", "avk", "sm_err_cov",
                            "n_err_cov"]:
                    results[key] = np.loadtxt(
                        self.workdir+"/"+key+".dat", skiprows=2)[:]
                for key in ["vcd", "dof", "flag", "counter"]:
                    results[key] = np.loadtxt(self.workdir+"/"+key+".dat")
                results["combierror"] = np.sqrt(
                    results["noiseerror"]**2 + results["smootherror"]**2 +
                    (results["profile"] * self.systerrfraction)**2)
                try:
                    scd = np.loadtxt(
                        self.workdir + '/SCD_table.dat', skiprows=2)[:, :4]
                except IndexError:
                    scd = np.loadtxt(
                        self.workdir + '/SCD_table.dat', skiprows=2)[:4]
                    scd = scd.reshape([1, scd.size])
                    TM.write_error("only one off-axis" +
                                   str(scd.shape), exc="", error=False)
                results["Sa_scaling"] = self.sa_scaling
                results["Sa_corrlen"] = self.sa_correlation_length
                results["avk_col"] = results["avk"].sum(axis=1)
                results["vcd_err"] = self.make_total_err(
                    results["sm_err_cov"], results["n_err_cov"],
                    results["vcd"])
            except (OSError, ValueError) as exc:  # mmf failed for some reason
                if TM.GSETTINGS["DEBUG"]:
                    if isinstance(exc, OSError):
                        TM.write_error("MMF failed tg in " + str(self.workdir),
                                       exc, False, False)
                    elif isinstance(exc, ValueError):
                        TM.write_error("MMF no tg attempted, ae nan in " +
                                       str(self.workdir), exc, False, False)
                results["profile"] = np.ones(self.totgrid_layers)*np.nan
                results["apriori"] = np.ones(self.totgrid_layers)*np.nan
                results["noiseerror"] = results["profile"]
                results["smootherror"] = results["profile"]
                results["combierror"] = results["profile"]
                scd = np.ones([self.numang, 4])*np.nan
                results["avk"] = np.ones(
                    [self.totgrid_layers, self.totgrid_layers])*np.nan
                results["sm_err_cov"] = np.copy(results["avk"])
                results["n_err_cov"] = np.copy(results["avk"])
                for key in results.oneD:
                    results[key] = np.nan
                results["vcd_err"] = np.nan
                results["avk_col"] = np.nan * np.ones(self.totgrid_layers)
        results["el_ang"] = scd[:, 0]
        results["scd_sim"] = scd[:, 1]
        results["scd_meas"] = scd[:, 2]
        results["scd_meas_err"] = scd[:, 3]
        results["scan_number"] = self.scan
        results["ssa"] = self.ssa
        results["asy"] = self.asy
        if np.isnan(results["vcd_err"]):
            results["flag"] = SETDICT["flags"]["fail"]
            results["rms_dscd"] = np.nan
            results["rms_r_dscd"] = np.nan
        else:
            results["flag"], results["rms_dscd"], results["rms_r_dscd"] = self.new_make_flag(
                results["flag"], results["scd_sim"], results["scd_meas"],
                results["scd_meas_err"], results["dof"], results["el_ang"])
        clean_mmf(self.workdir)
        del results["el_ang"]
        return results

    def make_total_err(self, sm_cov, noise_cov, vcd):
        '''
        Adding up the error matrices to construct total error on vcd

        This takes into account the error contributions from the noise error
        covariance matrix noise_cov, calculated as GSmG.T smoothing error
        covariance matrix sm_cov, calcualted as AkSaAk.T and a systematic error
        fraction of the retrieved vcd. This systematic error fration is set in
        the configuration fileIn theory, the smoothing and noise error
        covariance matrices have to be simply summed up in order to give the
        integrated error. However, in case of the smoothing matrix, given that
        the Sa is a constructed one and not representing the real variability
        in finer structure, this would result in a underestimation, hence, only
        the trace is used for the smoothing error contribution. This might
        change later.

        Parameters:
        ----------
        sm_cov: 2D array
            smoothing covariance matrix in partial column
        noise_cov: 2D array
            noise covariance matrix in partial column
        vcd: scalar
            vcd (or aod in case of aerosol) in the same units as sm_/noise_cov

        Returns:
        ---------
        totalerror: scalar
            the total vcd error
        '''
        noisecontribution = noise_cov.sum()
        smoothingcontribution = np.trace(sm_cov).sum()
        systematiccontr = (self.systerrfraction*vcd)**2
        totalerror = np.sqrt(
            noisecontribution + smoothingcontribution + systematiccontr)
        return totalerror

    @staticmethod
    def new_make_flag(flag, sim, meas, err, dof, elang):
        '''
        form a flag taking the rms, dof and convergence into account

        Parameters:
        -----------
        flag: integer
            0, 1 or 2 returned by fortran code. 0 ok, 1 warn, 2 not converged
        sim: 1D array
            simulated dscd
        meas: 1D array
            measured dscd
        err: 1D array
            measured dscd error
        dof: float
            degree of freedom
        elang: 1D aray
            levation angles in degree

        Returns:
        --------
        flag: int
            integer containing sum or errors as in SETDICT["flags"]
        rms: float
            rms of measured and modeled in measured error units
        rms_15: float
            as rms but only up to el ang SETDICT["settings"]["elang_limit"]
        '''
        if flag >= SETDICT["flag_limits"]["converged"]:
            flag = SETDICT["flags"]["converged"]
        elif flag >= SETDICT["flag_limits"]["converged_warn"]:
            flag = SETDICT["flags"]["converged_warn"]
        else:
            flag = 0
        if dof < SETDICT["flag_limits"]["dof"]:
            flag = flag + SETDICT["flags"]["dof"]
        elif dof < SETDICT["flag_limits"]["dof_warn"]:
            flag = flag + SETDICT["flags"]["dof_warn"]
        elang[np.isnan(elang)] = 999999
        # this was really just for testing
        if "elang_limit" not in SETDICT["settings"].keys():
            mask = elang < 100.
        else:
            mask = elang < SETDICT["settings"]["elang_limit"]
        err_min = np.fmax(err, SETDICT["settings"]["minimum_frac_error"]*meas)
        rms = np.sqrt(np.nanmean(((sim-meas)/err_min)**2))
        # only test purposes
        rms_15 = np.sqrt(np.nanmean(((sim[mask]-meas[mask])/err_min[mask])**2))
        # since mask is not masking anything if ["elang_limit"] not set
        # the below just returns rms
        min_rms = np.nanmin([rms, rms_15])
        if min_rms > SETDICT["flag_limits"]["rms"]:
            flag = flag + SETDICT["flags"]["rms"]
        elif min_rms > SETDICT["flag_limits"]["rms_warn"]:
            flag = flag + SETDICT["flags"]["rms_warn"]
        return flag, rms, rms_15


def ini_mmf(mpath):
    '''
    mpath: string
        path to mmf configuration file mmf_settings.yml containing:
        SIG_DIC, flags, flag_limits, settings
    '''
    global SETDICT
    with open("/".join([mpath, "mmf_settings.yml"])) as fid:
        SETDICT = yaml.load(fid, Loader=yaml.Loader)
    return

class MMF_MASTER(object):
    '''
    Master class to run MMF stand alone or within FRM4DOAS framework

    After making an MMF_MASTER object, call routine run
    '''

    def __init__(self, settings, mdict, ancil):  # settings
        '''
        Initialize MMF

        Parameters:
        -----------
        settings: dict
            settings dict with keys: mmftempdir, mmfconfigs, mmfexe_ae,
            mmfexe_tg, mmfinputs, multiprocs, procs
        mdict: dictionarry
            dictionarry with inputs, keys are: RAA, SAA, SZA, TAA, el_ang,
            scan_idx, doy, year, o4used, ae, tg. The latter two are dicts.
        ancil: dict
            ancil data dict such as temperature and aerosol properties
        '''
        self.settings = settings
        self.mdict = mdict
        self.ancil = ancil
        self.do_mmf = True
        ini_mmf(settings["mmfconfigs"])  # setting setdict
        try:
            self.extra_error = settings["extra_error"]
        except KeyError:
            self.extra_error = False

    def run(self):
        '''
        Run MMF retrieval, main function to call, wrapper around run_mmf

        Returns:
        --------
        mmf_results: dict
            results from MMF in form of a dictionarry, each value is a
            MMFobject
        erc: integer
            outgoing errorcode
        '''
        erc = 0
        if not self.settings["do_ae"]:  # for now, this needs to be done
            TM.write_error("Aerosol loading not yet implemented", "")
            self.mdict["do_mmf"] = False
            erc = erc + 1
        else:
            try:
                mmf_aerosols = {}
                for aer in list(self.mdict["ae"].keys()):
                    mmf_aerosols[aer], local_erc = self.run_mmf("aerosol", aer)
                if local_erc == 0:
                    if self.settings["do_tg"]:
                        try:
                            mmf_tracegases = {}
                            for gas in list(self.mdict["tg"].keys()):
                                mmf_tracegases[gas], _ = self.run_mmf(
                                    ["tracegas", mmf_aerosols[
                                        self.mdict["o4used"][gas]]], gas)
                        except Exception as exc:
                            TM.write_error(
                                "MMF TG failed, continue without MMF TG", exc)
                            # self.mdict["do_mmf"] = False
                            erc = erc + 1
                            mmf_tracegases = {}
                    else:
                        mmf_tracegases = {}
                erc = erc + local_erc
            except Exception as exc:
                TM.write_error("MMF AE failed, continue without MMF ", exc)
                self.do_mmf = False
                erc = erc + 1
                mmf_aerosols = None
                mmf_tracegases = None
            mmf_results = self.output(mmf_aerosols, mmf_tracegases)
            if erc > 0:
                raise ValueError("MMF failed")
            return mmf_results, erc

    def run_mmf(self, inputs, xgas):
        '''
        wrapper for running mmf retrieval either tg or ae

        prepare the call to mmf for each scan using aerosol properties or
        default values.

        Parameters:
        -----------
        inputs: list
            containing at least a string "aerosol" or "tracegas". If "tracegas"
            also needs to have aerosol results as second entry in list.
        xgas: string
            name of gas (or o4 window) used for retrieval

        Returns:
        --------
        results: list of 2-tuple
            each row is one scan. The first entry of each tuple is the scan
            number. The second is a list of results. The length and contents
            depend on whether aerosol or tracegas retrieval was performed.
        local_erc: 1D array of int
            array of local error codes from single scan

        Comments:
        ---------
        In the current setup, each retrieval is run twice: Once with ae apriori
        normal, once with modified a priori. The gas a priori is not changed,
        but run with the two different ae results. This is used as an
        additional flag, consistence flag.
        '''
        # extract the inputs
        if "aerosol" in inputs:
            mgas = self.mdict["ae"][xgas]
            is_tracegas = False
        elif "tracegas" in inputs:
            mgas = self.mdict["tg"][xgas]
            is_tracegas = True
        if self.settings["multiprocs"]:
            task_queue = multiprocessing.Manager().Queue()
            done_queue = multiprocessing.Manager().Queue()
            if TM.GSETTINGS["DEBUG"]:
                TM.dwrite("doing multip")
        else:
            if TM.GSETTINGS["DEBUG"]:
                TM.dwrite("doing single processing")
        windowname = (str(xgas)).upper()
        if "NO2" in windowname.upper():
            name = "NO2"
        elif "O4" in windowname.upper():
            name = "O4"
        elif "HCHO" in windowname.upper():
            name = "HCHO"
        elif "HONO" in windowname.upper():
            name = "HONO"
        results_all = {}
        runlist = ["NORMAL", "LOW_APRI"]
        for run in runlist:
            reslist = []
            mmf = []
            for ii, scan in enumerate(self.mdict["scan_idx"]):
                gasdict = {"dscd": mgas["dscd"][ii],
                           "dscd_err": mgas["dscd_err"][ii],
                           "wavelength": mgas["wavelength"],
                           "name": name}
                # This was only for testing
                # if self.extra_error:  # only if it was not included already
                #    gasdict["dscd_err"][2:] = np.sqrt(
                #        np.nanmax(mgas["dscd_extra_err"][ii][:2]
                #                  )**2+mgas["dscd_err"][ii][2:]**2)
                if is_tracegas:
                    # check for scan consistence, otherwise break out of loop?
                    if "LOW_APRI" in run:
                        if scan not in inputs[1]["LOW_APRI"]["scan_number"]:
                            print(scan, "is not performed twice")
                            continue
                        else:
                            idx = inputs[1]["LOW_APRI"]["scan_number"] == scan
                            ae_prof = inputs[1]["LOW_APRI"]["profile"][idx].flatten()
                            ae_tau = inputs[1]["LOW_APRI"]["aod"][idx]
                            if self.ancil["angstrom_exponent"][ii] != 0:
                                ae_tau = (
                                    (float(mgas["wavelength"]) / float(
                                        self.mdict["ae"][
                                            self.mdict["o4used"][xgas]]["wavelength"])
                                    )**(-self.ancil["angstrom_exponent"][ii])*ae_tau)
                    else:
                        if scan == inputs[1]["scan_number"][ii]:
                            if TM.GSETTINGS["DEBUG"]:
                                TM.dwrite("scan "+str(
                                    inputs[1]["scan_number"][ii])+ " is ok")
                        else:
                            raise ValueError("scan index inconsitence, tg: " +
                                             str(scan) + " aerosol: " + str(
                                                 inputs[1]["scan_number"][ii]))
                        ae_prof = inputs[1]["profile"][ii]
                        ae_tau = inputs[1]["aod"][ii]
                        if self.ancil["angstrom_exponent"][ii] != 0:
                            ae_tau = (
                                (float(mgas["wavelength"]) / float(
                                    self.mdict["ae"][
                                        self.mdict["o4used"][xgas]]["wavelength"])
                                )**(-self.ancil["angstrom_exponent"][ii])*ae_tau)
                else:
                    ae_prof = None
                    if self.ancil["aodapriori"] is None:
                        ae_tau = TM.DEFAULTS["aod_apriori"]
                    else:
                        if self.ancil["aodapriori"].ndim == 1:
                            ae_tau = self.ancil["aodapriori"]
                        else:
                            ae_tau = np.interp(
                                gasdict["wavelength"], self.ancil["aewav"],
                                self.ancil["aodapriori"][ii])
                    if "LOW_APRI" in run:
                        ae_tau = ae_tau*SETDICT["settings"]["test_apriori_aod_factor"]
                if (self.ancil["aeasym"]).ndim == 1:
                    asy = self.ancil["aeasym"][ii]
                else:
                    asy = np.interp(gasdict["wavelength"], self.ancil["aewav"],
                                    self.ancil["aeasym"][ii])
                if (self.ancil["aessa"]).ndim == 1:
                    ssa = self.ancil["aessa"][ii]
                else:
                    ssa = np.interp(gasdict["wavelength"], self.ancil["aewav"],
                                    self.ancil["aessa"][ii])
                tph = np.transpose(
                    np.stack([self.ancil["TPheight"], self.ancil["pressure"][ii],
                              self.ancil["temperature"][ii]]))
                mmf.append(MMF(
                    self.settings["mmftempdir"], tph.astype(str),
                    str(self.mdict["z_detector"]), scan,
                    self.settings["mmfinputs_specific"], gasdict,
                    self.mdict["RAA"][ii], self.mdict["SZA"][ii],
                    self.mdict["el_ang"][ii], ae_prof, ae_tau,
                    xgas, self.settings["mmfinputs"], self.settings["mmfexe_tg"],
                    self.settings["mmfexe_ae"], asy, ssa,
                    self.ancil["surface_albedo"]))
                if self.settings["multiprocs"]:
                    if is_tracegas:
                        task_queue.put((mmf[-1].run_tg, scan))
                    else:
                        task_queue.put((mmf[-1].run_ae, scan))
                else:
                    if is_tracegas:
                        reslist.append(mmf[-1].run_tg())
                    else:
                        reslist.append(mmf[-1].run_ae())
            if self.settings["multiprocs"]:
                for ii2 in range(self.settings["procs"]):
                    task_queue.put("STOP")
                proclist = []
                for ii2 in range(self.settings["procs"]):
                    proclist.append(Process(target=processwrapper,
                                            args=(task_queue, done_queue)))
                    if TM.GSETTINGS["DEBUG"]:
                        TM.dwrite("started process "+str(ii2))
                    proclist[-1].start()
                for proc in proclist:
                    proc.join()
                    if TM.GSETTINGS["DEBUG"]:
                        TM.dwrite("joined " + str(proc))
                for ii2 in range(done_queue.qsize()):
                    mres = done_queue.get()
                    reslist.append(mres)
            mmf_result = MMFresults()
            local_erc = mmf_result.join(reslist)
            results_all[run] = mmf_result
        if not is_tracegas:
            extra = self.check_consistence(
                results_all["NORMAL"], results_all["LOW_APRI"])
        else:
            if len(results_all["LOW_APRI"]) > 0:
                extra = results_all["LOW_APRI"]
            else:
                extra = {"aod": [], "profile": [], "scan_number": []}
        mmf_result = results_all["NORMAL"]
        mmf_result["LOW_APRI"] = extra
        return mmf_result, local_erc

    @staticmethod
    def check_consistence(normal, test):
        '''
        Check consistence btw ae runs normal and test to decide if tg run twice

        Parameters:
        -----------
        normal: dict
            containing at least keys: aod, profile, scan_number of orig ae run
        test: dict
            as normal, but with the varied a priori settings

        Returns:
        --------
        test: dict
            as input test, but only the ones that are not consistent
        '''
        fraction = np.maximum(
            normal["aod"]/test["aod"], test["aod"]/normal["aod"])
        mask = fraction > SETDICT["flag_limits"]["stability_aod_warn"]
        if mask.sum() > 0:
            test = {key: test[key][mask] for key in test.keys() if
                    key in ["aod", "profile", "scan_number"]}
        else:
            test = {"aod": [], "profile": [], "scan_number": []}
        return test

    def output(self, mmf_aerosols, mmf_tracegases):
        '''
        Collect data into a single results dictionarry with correct keys

        Parameters:
        -----------
        mmf_aerosols: dictionarry
            dictionarry holding aerosol results
        mmf_tracegases: dictionarry
            dictionarry holding tracegas results

        Returns:
        -------
        mmfout: dictionarry
            results dictionarry with keys as nc paths and values as dict
            containing the value, FillValue, units, description, name
        '''
        with open(os.path.join(self.settings["mmfconfigs"],
                               SETDICT["settings"]["outputconfigname"]
                              )) as fid:
            if float(yaml.__version__.split(".")[0]) < 5:
                mmfdict = yaml.load(fid)
            else:
                mmfdict = yaml.load(fid, Loader=yaml.FullLoader)
        mmfout = {}
        for aer in mmf_aerosols.keys():
            aecheck = mmf_aerosols[aer].pop("LOW_APRI")
            scan_number = mmf_aerosols[aer].pop("scan_number")
            for ii1, number in enumerate(scan_number):
                if number in aecheck["scan_number"]:
                    idx = aecheck["scan_number"] == number
                    test = aecheck["aod"][idx]
                    normal = mmf_aerosols[aer]["aod"][ii1]
                    masking = np.maximum(test/normal, normal/test)
                    if masking > SETDICT["flag_limits"]["stability_aod"]:
                        mmf_aerosols[aer]["flag"][ii1] = (
                            mmf_aerosols[aer]["flag"][ii1]
                            + SETDICT["flags"]["stability"])
                    elif masking > SETDICT["flag_limits"]["stability_aod_warn"]:
                        mmf_aerosols[aer]["flag"][ii1] = (
                            mmf_aerosols[aer]["flag"][ii1] +
                            SETDICT["flags"]["stability_warn"])
            # change the current flag to flag_detail and
            # form the flag from the flag_detail
            mmf_aerosols[aer]["flag_detail"] = np.copy(
                mmf_aerosols[aer]["flag"][:])
            mmf_aerosols[aer]["flag"][(mmf_aerosols[aer]["flag"] > 0) &
                                      (mmf_aerosols[aer]["flag"] <
                                       SETDICT["settings"]["is_error"])] = 1
            mmf_aerosols[aer]["flag"][
                mmf_aerosols[aer]["flag"] >= SETDICT["settings"]["is_error"]
                ] = 2
            orig_layerthickness = mmf_aerosols[aer].pop("layerthickness")
            version = mmf_aerosols[aer].pop("version")
            mmfheights = mmf_aerosols[aer].pop("middleheights")
            if "height_boundaries_m" in self.mdict:
                desired_heights = self.mdict["height_boundaries_m"]
            else:
                desired_heights = np.array(
                    TM.REALDIMS["HEIGHTGRIDKM"])*TM.CONV["KM2M"]
            layerthickness = np.abs(np.diff(desired_heights))
            if desired_heights[0] > desired_heights[-1]:
                desired_heights = desired_heights[1:] + layerthickness * 0.5
            else:
                desired_heights = desired_heights[:-1] + layerthickness * 0.5
            desired_heights = desired_heights + (
                self.mdict["z_detector"]*TM.CONV["KM2M"])
            for key in mmf_aerosols[aer]:
                # consult the aerosol key and find correct name and stuff
                try:
                    mentry = copy.deepcopy(mmfdict[key+"_aer"])
                except KeyError:
                    # this means that it is something that should not be saved
                    # in the output, so just jump to the next output
                    continue
                mkey = mentry["name"].replace("aer/", aer.upper() + "/")
                mdata = mmf_aerosols[aer][key]
                if "DIM_LAYERS" in mentry["dimensions"]:
                    if "aircol" in key:
                        mdata = mdata / orig_layerthickness
                    mdata = interpolate_data_right(
                        mentry, mdata, mmfheights, desired_heights)
                    # TODO: I should probably use log for pressure and aircol
                    if "aircol" in key:
                        mdata = mdata * layerthickness
                mentry["value"] = mdata
                mentry["name"] = mentry["name"].replace("aer/", "")
                mmfout[mkey] = mentry
            for mkey in ["desired_heights", "layerthickness"]:
                if mkey in mmfdict.keys():
                    mentry = copy.deepcopy(mmfdict[mkey])
                    mentry["value"] = locals()[mkey] # in m
                    mmfout[mentry["name"]] = mentry
            if "FROM_ORIG" in mmfdict.keys():
                for mkey in mmfdict["FROM_ORIG"]:
                    mentry = copy.deepcopy(mmfdict["FROM_ORIG"][mkey])
                    mentry["value"] = self.mdict[mkey] # in m
                    mmfout[mentry["name"]] = mentry
        if self.settings["do_tg"]:
            for tg in mmf_tracegases.keys():
                scan_number = mmf_tracegases[tg].pop("scan_number")
                tracegastest = mmf_tracegases[tg].pop("LOW_APRI")
                for ii1, number in enumerate(scan_number):
                    if number in tracegastest["scan_number"]:
                        idx = tracegastest["scan_number"] == number
                        test = tracegastest["vcd"][idx]
                        normal = mmf_tracegases[tg]["vcd"][ii1]
                        masking = np.maximum(test/normal, normal/test)
                        if masking > SETDICT["flag_limits"]["stability_vcd"]:
                            mmf_tracegases[tg]["flag"][ii1] = (
                                mmf_tracegases[tg]["flag"][ii1] +
                                SETDICT["flags"]["stability"])
                        elif masking > SETDICT["flag_limits"
                                              ]["stability_vcd_warn"]:
                            mmf_tracegases[tg]["flag"][ii1] = (
                                mmf_tracegases[tg]["flag"][ii1] +
                                SETDICT["flags"]["stability_warn"])
                # change the current flag to flag_detail and
                # form the flag from the flag_detail
                mmf_tracegases[tg]["flag_detail"] = np.copy(
                    mmf_tracegases[tg]["flag"][:])
                mmf_tracegases[tg]["flag"][(mmf_tracegases[tg]["flag"] > 0) & (
                    mmf_tracegases[tg]["flag"] <
                    SETDICT["settings"]["is_error"])] = 1
                mmf_tracegases[tg]["flag"][
                    (mmf_aerosols[self.mdict["o4used"][tg]]["flag"] > 1) &
                    (mmf_tracegases[tg]["flag"] < SETDICT["settings"]["is_error"])
                                          ] = 1
                mmf_tracegases[tg]["flag"][
                    mmf_tracegases[tg]["flag"] >= SETDICT["settings"]["is_error"]
                    ] = 2
                for key in mmf_tracegases[tg]:
                    try:
                        mentry = copy.deepcopy(mmfdict[key+"_tg"])
                    except KeyError:
                        # means that it is something that should not be saved
                        # in the output, so just jump to the next output
                        continue
                    # consult the tg key and find the correct name and stuff
                    mkey = mentry["name"].replace("tg/", tg.upper() + "/")
                    mdata = mmf_tracegases[tg][key]
                    if "DIM_LAYERS" in mentry["dimensions"]:
                        mdata = interpolate_data_right(
                            mentry, mdata, mmfheights, desired_heights)
                    mentry["value"] = mdata
                    mentry["name"] = mentry["name"].replace("tg/", "")
                    mmfout[mkey] = mentry
        mmfout["attributes"] = {}
        mmfout["attributes"]["version"] = version
        return mmfout

def interpolate_data_right(mentry, mdata, mmfheights, desired_heights):
    '''
    perform data interpolation of mdata on mmfheights to desired_heights

    Parameters:
    -----------
    mentry: dict
        needs to have key "dimensions", these are dimensions of mdata
    mdata: array like
        input data to be interpolated
    mmfheights: array like
        array of heights on which the input data is
    desired_heights: array like
        array of heights to interpolate on

    Returns:
    --------
    mdata: array like
        input data interpolated on desired_heights, if there was a dimension
        that corresponds to a height dimension.
    '''
    entrs = np.where(np.array(mentry["dimensions"]) == "DIM_LAYERS")[0]
    mentrydimlen = len(mentry["dimensions"])
    # since calling interp2d returnes values on ascending grid,
    # eventhough it is actually not ascending, need to re-do the order
    order = np.argsort(desired_heights)
    o_order = np.argsort(order)
    if mentrydimlen == 1:
        mdata = np.interp(desired_heights, mmfheights[::-1], mdata[::-1])
    if mentrydimlen == 2:
        if len(entrs) == 1:
            if entrs == 0:
                mdn = []
                for col in mdata.T:
                    mdn.append(np.interp(
                        desired_heights, mmfheights[::-1], col[::-1]))
                mdata = (np.array(mdn)).T
            elif entrs == 1:
                mdn = []
                for row in mdata:
                    mdn.append(np.interp(
                        desired_heights, mmfheights[::-1], row[::-1]))
                mdata = np.array(mdn)
        else:
            spf = interp2d(mmfheights, mmfheights, mdata)
            mdata = spf(desired_heights, desired_heights)
            # need to re-order
            mdata = mdata[o_order, :]
            mdata = mdata[:, o_order]
            # 2D interpolation
    if mentrydimlen == 3:
        if len(entrs) == 1:
            raise ValueError(
                "dimension of supplied data =3, but dimlayer only appear once")
            # 1D interpolation, but for each row, column or height
        if len(entrs) == 2:
            mdn = []
            if 1 in entrs and 2 in entrs:
                for mslice in mdata:
                    spf = interp2d(mmfheights, mmfheights, mslice)
                    interm = spf(desired_heights, desired_heights)
                    interm = interm[o_order, :]
                    interm = interm[:, o_order]
                    mdn.append(interm)
                mdata = np.array(mdn)
            elif 0 in entrs and 1 in entrs:
                for mslice in mdata.T:
                    spf = interp2d(mmfheights, mmfheights, mslice)
                    interm = spf(desired_heights, desired_heights)
                    interm = interm[o_order, :]
                    interm = interm[:, o_order]
                    mdn.append(interm)
                mdata = (np.array(mdn)).T
            else:
                raise ValueError(
                    "dimension of supplied data =3, layerdim at 0 and 2")
    return mdata

def mmf_stand_alone(settings, ancil, mdict, outputname):
    '''
    Routine to call if mmf is run as a stand alone code

    Parameters:
    -----------
    settings: dict
        settings dictionarry with keys: do_ae, do_tg, lower_el_ang_limit
        mmfconfigs, mmfexe_ae, mmfexe_tg, mmfinputs, mmfinputs_specific,
        mmftempdir, configuration, multiprocs, procs, dimension_names
    ancil: dict
        ancilary data dictionarry with keys: TPheight, aeasym, aessa,
        aewav, angstrom_exponent, aodapriori, pressure, surface_albedo,
        temperature
    mdict: dict
        measurement data dictionarry with keys: RAA, SZA, RAA, (TAA), ae, tg
        (doy), el_ang, (year), o4used, scan_idx, z_detector
    outputnames: str
        name of output netCDF4 file to be created
    '''
    TM.initialize(settings["configuration"])
    mmf = MMF_MASTER(settings, mdict, ancil)
    mmf, erc = mmf.run()
    mmf_new = {}
    for key in mmf:
        try:
            mmf_new[key] = TM.translate_var(mmf[key], settings["dimension_names"])
        except:
            mmf_new[key] = mmf[key]
    dims = {
        settings["dimension_names"]["DIM_LAYERS"]: TM.REALDIMS["NUM_GRID"],
        settings["dimension_names"]["DIM_SCAN_NAME"]: len(mdict["scan_idx"]),
        settings["dimension_names"]["DIM_ANGLE_NAME"]: TM.REALDIMS["MAXAXIS"]}
    mmf_exp = yn.expand_keys(mmf_new)
    yn.dict_to_nc(mmf_exp, outputname, dims)
    return erc

def test_run(settings_loc, ancil_loc, mdict_loc):
    '''
    Test the program with the supplied test files.

    Parameters:
    -----------
    settings_loc: str
        path to settings yaml file
    ancil_loc: str
        path to ancil yaml file
    mdict_loc: str
        path to measurement yaml file
    '''
    start = time.time()
    with open(settings_loc) as fid:
        settings = yaml.load(fid, Loader=yaml.Loader)
    with open(ancil_loc) as fid:
        ancil = yaml.load(fid, Loader=yaml.Loader)
        for key in ancil:
            if hasattr(ancil[key], '__len__'):
                ancil[key] = np.array(ancil[key]).astype(float)
    with open(mdict_loc) as fid:
        mdict = yaml.load(fid, Loader=yaml.Loader)
        for key in mdict:
            if isinstance(mdict[key], dict):
                for ky1 in mdict[key]:
                    try:
                        for ky2 in mdict[key][ky1]:
                            mdict[key][ky1][ky2] = np.array(
                                mdict[key][ky1][ky2]).astype(float)
                    except:
                        pass
            else:
                try:
                    mdict[key] = np.array(mdict[key])
                except:
                    pass
    ecr = mmf_stand_alone(settings, ancil, mdict, "test.nc")
    print("error code", ecr)
    stop = time.time()
    print(stop-start)

if __name__ == "__main__":
    if len(sys.argv) == 4:
        (SETTINGS_F, ANCIL_F, MDICT_F) = sys.argv[1:4]
    else:
        SETTINGS_F = "TESTDATA/SETTINGS.yml"
        ANCIL_F = "TESTDATA/ANCIL.yml"
        MDICT_F = "TESTDATA/MEAS.yml"
    test_run(SETTINGS_F, ANCIL_F, MDICT_F)
