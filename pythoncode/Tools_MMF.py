'''
Collection of methods and classes to use with MMF and more generally FRM4DOAS

Author: Martina M. Friedrich

Date: Nov 2019
'''
import inspect
import sys
import traceback
import datetime
import re
import yaml
import netCDF4
import numpy.ma as ma
import numpy as np
from cf_units import Unit

DIMENSIONS = {}
FILLVALS = {}
GSETTINGS = {}
REALDIMS = {}
ERRCODES = {}
DEFAULTS = {}
SYSERRFRAC = {}

def initialize(configname):
    '''
    function to initialize settings via a configuration file

    Parameters:
    -----------
    configname: str
        path to configuration file name in yaml format with keys:
            FILLVALS, GSETTINGS, REALDIMS, ERRCODES, DEFAULTS, SYSERRFRAC
    '''
    global DIMENSIONS, FILLVALS, GSETTINGS, REALDIMS, \
        ERRCODES, DEFAULTS, SYSERRFRAC
    if isinstance(configname, str):
        fid = open(configname)
        if float((yaml.__version__).split(".")[0]) < 5:
            conf = yaml.load(fid)
        else:
            conf = yaml.load(fid, Loader=yaml.FullLoader)
    else:
        conf = configname
    FILLVALS = conf["FILLVALS"]
    GSETTINGS = conf["GSETTINGS"]
    REALDIMS = conf["REALDIMS"]
    ERRCODES = conf["ERRCODES"]
    DEFAULTS = conf["DEFAULTS"]
    SYSERRFRAC = conf["SYSERRFRAC"]
    return


CONV = {
    "M2KM": 0.001,
    "DAY2SEC": 24*360.,
    "KM2M": 1000,
    "M2CM": 100,
    "Pa2hPa": 0.01,
    "KM2CM": 100000.,
    "hPa2Pa": 100.}

CONST = {
    "k_B": 1.380649e-23,
    "g_grav": 9.80665,  # m/s**2
    "R_d": 287.  # m**2/s**2/K  specific gas constant R/Ma
    }

def translate_var(idn_name, dimdic):
    '''
    translate the dimension and fillvalue names for a variable

    Parameters:
    -----------
    idn_name: dict
        yaml_netcdf4 dictionarry for a variable
    dimdic: dict
        dict containing key "dimensions" with keys dimension names

    Returns:
    --------
    idn_name: dict
        yaml_netcdf4 dictionarry for variable with translated dims and fillvals
    '''
    dim = idn_name["dimensions"]
    idn_name["dimensions"] = tuple([dimdic[entr]for entr in dim if entr != ""])
    try:
        fillval = idn_name["_FillValue"]
        idn_name["_FillValue"] = FILLVALS[fillval]
    except KeyError:
        try:
            fillval = idn_name["attributes"]["_FillValue"]
            idn_name["_FillValue"] = FILLVALS[fillval]
        except KeyError:
            pass
    return idn_name

def set_attributes(handle, attrs):
    '''
    set attributes in attrs on a netCDF4 handle

    Parameters:
    -----------
    handle: netCDF4 group/ file or variable handle
        handle to add attributes to
    attrs: dict
        dictionarry of attributes to add
    '''
    handle.setncatts({key: attrs[key] for key in attrs.keys()})
    return

def get_attributes(handle):
    '''
    get all attributes attached to the file_handle of nc file

    Parameters:
    -----------
    file_handle: netCDF4 group/ file/ variable handle
        handle to get the attributes from

    Returns:
    --------
    mdict: dict
        dictionarry of attributes that were attached to the handle
    '''
    keys = handle.ncattrs()
    mdict = {}
    for key in keys:
        mdict[key] = handle.getncattr(key)
    return mdict

def get_group_variables(group_handle):
    '''
    get content of a group from a nc group or file handle

    Parameters:
    -----------
    group_handle: netCDF4 group/ file handle
        any group/ file handle

    Returns:
    -------
    allvars: dict
        dict of handle variables
    allgroups: dict
        dict of handle groups
    '''
    all_vars = group_handle.variables
    all_groups = group_handle.groups
    allvars = {}
    allgroups = {}
    for var in all_vars:
        allvars[var] = group_handle.variables[var]
    for group in all_groups:
        allgroups[group] = group_handle.groups[group]
    return allvars, allgroups

def copy_whole_group(group_old, group_new, limit=None):
    '''
    recursive function to copy a whole nc group

    Parameters:
    -----------
    group_old: group handle/ file handle
        of old group
    group_new: group handle/ file handle
        of new group
    [limit: integer
        this is only to copy qdoas files, to cut off scans for tests]
    '''
    for dname, the_dim in group_old.dimensions.items():
        if not the_dim.isunlimited():
            if limit is None or dname not in "n_alongtrack":
                group_new.createDimension(dname, len(the_dim))
            else:
                group_new.createDimension(dname, limit)
        else:
            group_new.createDimension(dname, None)
    for key in group_old.groups.keys():
        group_new.createGroup(key)
        copy_whole_group(group_old[key], group_new[key], limit=limit)
    set_attributes(group_new, get_attributes(group_old))
    for vkey, vin in group_old.variables.items():
        try:
            outvar = group_new.createVariable(
                vkey, vin.datatype, vin.dimensions)
        except ValueError:
            for dim in vin.dimensions:
                mff = group_old.parent
                while True:
                    if dim in mff.dimensions:
                        group_new.createDimension(dim, len(mff.dimensions[dim]))
                        break
                    else:
                        mff = mff.parent
            outvar = group_new.createVariable(
                vkey, vin.datatype, vin.dimensions)
        for k in vin.ncattrs():
            try:
                outvar.setncatts({k: vin.getncattr(k)})
            except:
                outvar.setncattr_string(k, vin.getncattr(k))
        if limit is None or ("n_alongtrack" not in vin.dimensions):
            try:
                outvar[:] = vin[:]
            except:
                print("failed to copy ", vkey)
        else:
            outvar[:] = vin[:limit]
    return

def reorder_field(index_field, field, olen, fillval, mdtype):
    '''
    use an index field (of shape new field with indices of old) to reorder

    Parameters:
    ----------
    index_field: 2D array of integers
        shape of new data, entries are indices of old format
    field: 1D or 2D array
        old data to put in the new format
    olen: integer
        length of original 1D data, to check if field is transposed
    fillval: float/ integer
        FillValue to use
    mdtype: string
        representing dtype

    Returns:
    --------
    new_fields: 2D or 3D array
        data input field in new order: (orig extra dim x) scans x MAXAXIS
    '''
    if field.ndim > 1:
        if field.shape[0] == olen:
            field = np.transpose(field)
            was_t = True
        else:
            was_t = False
    else:
        field = [field]
        was_t = False
    new_fields = []
    for mfield in field:
        new_field = np.full(index_field.shape, fillval, mdtype)

        for i2i, line in enumerate(index_field):
            for j2j, entry in enumerate(line):
                if entry > -1:  # was !=-1
                    new_field[i2i, j2j] = mfield[index_field[i2i, j2j]]
        # for the moment, some of the int fields in QDOAS out, have -1
        # where they should have the fillval. Fix this here.
        if mdtype == "int32":
            new_field[new_field == -1] = fillval
        new_fields.append(new_field)
    new_fields = np.array(new_fields)
    if was_t:
        new_fields = np.moveaxis(new_fields, 0, -1)
    return np.squeeze(new_fields)

def dwrite(message):
    '''
    test function to print

    Parameters:
    -----------
    message: str
        message to print now
    '''
    print(" " + message, flush=True)
    return

def write_error(message, exc="", error=True, printtrace=True):
    '''
    method to print error/ warning message in desired form + produce traceback

    [file.py: function] ERROR: message exc (The error raised)
    ---------------------------------------------------------
    TRACE:
    traceback
    ---------------------------------------------------------

    Parameters:
    -----------
    message: string
        manual details to print in error, main message
    exc: string/ Error object
        last produced error message
    [error: logical]
        if True, report error if False, report warning, default True
    [printtrace: logical]
        if True print error trace otherwise not print error trace, default True
    '''
    try:
        frame = inspect.getouterframes(inspect.currentframe(), 2)[1]
        where = ": ".join([frame.filename, frame.function])
    except:
        where = sys._getframe().f_back.f_code.co_filename
    if error:
        print("[", where, "] ERROR: " + message + " " + str(exc), flush=True)
    else:
        print("[", where, "] WARNING: " + message + " " + str(exc), flush=True)
    if printtrace:
        if sys.exc_info()[0]:
            print("_________________________________________________________")
            print("TRACE:", flush=True)
            print(traceback.format_exc(), flush=True)
            print("_________________________________________________________")
    return

def versionconvert(version):
    '''
    convert the "free" version into a geoms valid version. Maybe not needed?

    Paramters:
    ----------
    version: str, float or int
        version name

    Returns:
    --------
    version: 4 digit string
        only numerical digits

    This only keeps version numbers until the 2nd decimal.
    '''
    if isinstance(version, str):
        try:
            version = float(re.sub("[^0-9.]", "", version))
        except ValueError:
            version = str(re.sub("[^0-9]", "", version))
            return version
    version = str(int(version*100)).zfill(4)
    return version

def fracday2int(day, year):
    '''
    Transforms fractional day and year into integer YYYYMMDDhhmm

    Parameters:
    -----------
    day: float
        fractional day of year
    year: int
        year

    Returns:
    --------
    dateint: integer
        integer of the form YYYYMMDDhhmm
    '''
    mdat = datetime.datetime(
        year, 1, 1, 0, 0, 0) + datetime.timedelta(day - 1.0)
    yrs = str(mdat.year).zfill(4)
    month = str(mdat.month).zfill(2)
    day = str(mdat.day).zfill(2)
    hour = str(mdat.hour).zfill(2)
    minute = str(mdat.minute).zfill(2)
    _ = str(mdat.second).zfill(2)
    dateint = int("".join([yrs, month, day, hour, minute]))
    return dateint

def local_solar_time(dto, longit):
    '''
    routine to calculate local solar time

    This follows documentation on
    [www.esrl.noaa.gov](http://www.esrl.noaa.gov/gmd/grad/solcalc/solareqns.PDF)

    Parameters:
    ------------
    dto: datetime.datetime object
        UT time
    longit: float
        longitude to calculate local solar time at

    Returns:
    --------
    solar_time: datetime.datetime object
        local solar time
    '''
    if hasattr(dto, "__len__"):
        dto = np.asarray(dto)
        oshape = dto.shape
        dto = dto.flatten()
    else:
        oshape = ()
        dto = [dto]
    solar_time = []
    for ddt in dto:
        gamma = 2.0*np.pi/365*(ddt.timetuple().tm_yday-1+float(ddt.hour-12)/24)
        eqtime = 229.18 * (
            0.000075 +
            0.00186800 * np.cos(1.*gamma) -
            0.032077 * np.sin(1.*gamma) -
            0.014615 * np.cos(2.*gamma) -
            0.040849 * np.sin(2.*gamma))
        time_offset = eqtime + 4. * longit  # this is in minutes
        tst = (ddt.hour * 60 + ddt.minute + ddt.second / 60 + time_offset)
        solar_time.append(datetime.datetime.combine(
            ddt.date(), datetime.time(0)) + datetime.timedelta(minutes=tst))
    solar_time = np.asarray(solar_time)
    if oshape != ():
        solar_time = solar_time.reshape(oshape)
    return solar_time

def fracday2datetime(fracday, year):
    '''
    Transforms fractional day and year into datetime object

    Parameters:
    -----------
    fracday: float or 1D array
        fractional day of year, starting at 1 on 1.1. of year
    year: int or 1D array
        year

    Returns:
    --------
    datetime object
    '''
    '''
    Transforms fractional day and year into string YYYYMMDDThhmmssZ

    Parameters:
    -----------
    fracday: float or 1D array
        fractional day of year, starting at 1 on 1.1. of year
    year: int or 1D array
        year

    Returns:
    --------
    string: YYYYMMDDThhmmssZ
    '''

    if np.ndim(fracday) == 0 and np.ndim(year) == 0:
        mdat = datetime.datetime(
            year, 1, 1, 0, 0, 0) + datetime.timedelta(fracday - 1.0)
    else:
        mdat = []
        for day, yyr in zip(fracday, year):
            mdat.append(fracday2datetime(day, yyr))
    return mdat

def fracday2datestr(fracday, year):
    '''
    Transforms fractional day and year into string YYYYMMDDThhmmssZ

    Parameters:
    -----------
    fracday: float or 1D array
        fractional day of year, starting at 1 on 1.1. of year
    year: int or 1D array
        year

    Returns:
    --------
    string: YYYYMMDDThhmmssZ
    '''

    if np.ndim(fracday) == 0 and np.ndim(year) == 0:
        mdat = datetime.datetime(
            year, 1, 1, 0, 0, 0) + datetime.timedelta(fracday - 1.0)
        mdate = datetime2str(mdat)
    else:
        mdate = []
        for day, yyr in zip(fracday, year):
            mdate.append(fracday2datestr(day, yyr))
        if isinstance(year, ma.core.MaskedArray
                     ) or isinstance(year, ma.core.MaskedArray):
            mdate = ma.array(mdate)
        elif isinstance(year, np.ndarray) or isinstance(fracday, np.ndarray):
            mdate = np.array(mdate)
    return mdate

def days2sec(fracday):
    '''
    Transforms a fractional day delta into seconds

    Parameters:
    -----------
    fracday: float or array of float
        fractional days timedelta

    Returns:
    --------
    seconds: float or array of float
        timedelta in seconds
    '''
    seconds = fracday * 24. * 60. * 60.
    return seconds

def datetime2str(mdat):
    '''
    transforms datetime object to string of the form YYYYMMDDThhmmssZ

    Parameters:
    -----------
    mdat: datetime.datetime object

    Returns:
    --------
    datestr: string
        datestring of the form YYYYMMDDThhmmssZ
    '''
    year = str(mdat.year).zfill(4)
    month = str(mdat.month).zfill(2)
    day = str(mdat.day).zfill(2)
    hour = str(mdat.hour).zfill(2)
    minute = str(mdat.minute).zfill(2)
    secs = str(mdat.second).zfill(2)
    datestr = "".join([year, month, day, "T", hour, minute, secs, "Z"])
    return datestr

def datestr2list(datestring):
    '''
    transforms data str like 20180526T030445Z to tuple(YY, MM, DD, hh, mm, ss)

    Parameters:
    -----------
    datestring: sring
        date string of the form YYYYMMDDThhmmssZ

    Returns:
    --------
    tuple of 6 integers year, month, day, hour, minute, second
    '''
    year = int(datestring[:4])
    month = int(datestring[4:6])
    day = int(datestring[6:8])
    hour = int(datestring[9:11])
    minute = int(datestring[11:13])
    second = int(datestring[13:15])
    mdate = (year, month, day, hour, minute, second)
    return mdate

def datestr2datetime(datestring):
    '''
    transforms data string of form 20180526T030445Z to datetime object

    Parameters:
    -----------
    datestring: sring
        date string of the form YYYYMMDDThhmmssZ

    Returns:
    --------
    datetime object
    '''
    return datetime.datetime(*datestr2list(datestring))

def datetime2mjd2000(*mdatetime):
    '''
    convert datetime, datestr or tuple to seconds since January 1st 2000 0:00

    Parameters:
    -----------
    mdatetime: datetime.datetime object, datestring YYYYMMDDThhmmssZ or 2-tuple
        in case of 2-tuple: (fractional day, year)

    Returns:
    --------
    int: time since January 1st 2000, 0:00 in fractional days
    '''
    if isinstance(mdatetime[0], datetime.datetime):
        pass
    elif isinstance(mdatetime[0], str):
        mdatetime = datestr2datetime(mdatetime)
    elif np.isscalar(mdatetime[0]) and len(mdatetime) == 2 and np.isscalar(mdatetime[1]):
        try:
            mdatetime = fracday2datestr(mdatetime[0], mdatetime[1])
        except np.ma.core.MaskError:
            return np.nan
        mdatetime = datestr2datetime(mdatetime)
    else:
        if ma.is_masked(mdatetime[0]):
            return np.nan
        else:
            raise ValueError(
                "wrong input to datetime2mjd200: " +
                "no datetime, date string or (fracday, year)")
    try:
        mdate = (mdatetime-datetime.datetime(2000, 1, 1, 0, 0, 0))
        fracday = mdate.days + mdate.seconds/datetime.timedelta(1).total_seconds()
    except ma.core.MaskError:
        fracday = np.nan
    return fracday

def datetime2fracday(mdatetime):
    '''
    return the fractional day of year given a datetime object

    Parameters:
    -----------
    mdatetime: datetime object
        datetime to convert

    Returns:
    --------
    fday: float
        the fractional day of year startin at 1 on January 1st
    '''
    if hasattr(mdatetime, "__len__"):
        fday = []
        for mdd in mdatetime:
            fday.append(datetime2fracday(mdd))
        fday = np.asarray(fday)
    else:
        timedelta = mdatetime - datetime.datetime(mdatetime.year, 1, 1)
        fday = timedelta.days+timedelta.seconds/(24*60*60)+1
    return fday

def ymlload(myfile):
    with open(myfile) as fid:
        if float(yaml.__version__.split(".")[0]) >= 5:
            mdict = yaml.load(fid, Loader=yaml.Loader)
        else:
            mdict = yaml.load(fid)
    return mdict

def convert2mjd2000(doy, year):
    '''
    Funtion to convert fractional DOY and YEAR to Fractional day since 2000

    This function uses MJD2000 from doas_lib

    Parameters:
    -----------
    doy: float or 1D array of float
        day of year
    year: int or 1D array of int

    Returns:
    --------
    mydate: flaot or 1D array of float
        the number of fractional days since 01/01/2000

    '''
    if np.array(year).shape == () and np.array(doy).shape == ():
        mydate = datetime2mjd2000(doy, year)
    else:
        mydate = []
        if np.array(year).size <= 1:
            year = (np.ones(doy.shape)*year).astype(int)
        for yyy, ddd in zip(year, doy):
            try:
                temp = convert2mjd2000(ddd, yyy)
            except np.ma.core.MaskError:
                temp = np.nan
            mydate.append(temp)
        mydate = np.array(mydate)
    return mydate

def mjd2000_2_datestr(idate):
    '''
    converts seconds since Jan, 1st 2000 00:00:00 to datestr YYYYMMDDThhmmssZ

    Parameters:
    -----------
    idate: float
        time in fractional days since Januar 1st, 2000 at 00:00:00

    Returns:
    --------
    datestr: str
        date as YYYYMMDDThhmmssZ

    '''
    try:
        idate = datetime2str(
            datetime.datetime(2000, 1, 1)+datetime.timedelta(idate))
    except ValueError:
        idate = ""
    return idate

def unitcheck(mvar, wunit):
    '''
    Convert the variable to the desired Unit

    Parameters:
    -----------
    mvar: netCDF4 var or dict with key "units" or "attributes" and then "units"
        needs to have "units" attribute
    wunit: string or Unit object
        needs to represent a valid CF Unit

    Returns:
    -------
    mval: array-type
        value of netCDF4 variable in unit wunit
    '''
    if isinstance(mvar, netCDF4._netCDF4.Variable):
        mval = Unit(mvar.getncattr("units")).convert(mvar[:], wunit)
    elif isinstance(mvar, dict):
        if "units" in mvar.keys():
            mval = Unit(mvar["units"]).convert(mvar["value"], wunit)
        elif "attributes" in mvar.keys():
            if "units" in mvar["attributes"].keys():
                mval = Unit(mvar["attributes"]["units"]).convert(
                    mvar["value"], Unit(wunit))
            else:
                raise ValueError("no units found in attributes for unitcheck")
        else:
            raise ValueError("no units or attributes found for unitcheck")
    else:
        raise ValueError("unexpected type for unitcheck: " + str(type(mvar)))
    return mval

def get_erc(errocodes):
    '''
    get a list of the errocodes present in the sum of errorcodes

    currently, the error-code 256 for geoms can be present several times.
    Take care of this

    Parameters:
    -----------
    errcodes: integer
        sum of error codes

    Returns:
    -------
    eclist: list of integers
        list of error codes contained
    '''
    eclist = []
    test = np.copy(errocodes)
    for ectest in ([2**entry for entry in range(12)])[::-1]:
        if test - ectest >= 0:
            eclist.append(ectest)
            test = test - ectest
        else:
            continue
    return eclist

def erc_in_ercs(errocodes, specific):
    '''
    check whether an error code is in a sum of error codes

    Parameters:
    -----------
    errocodes: integer
        sum of potenz of 2 representing a cummulative error code
    specific: integer
        potenz of 2, representing a specific error

    Returns:
    --------
    isin: logical
        True if specific error code forms part of the error, False otherwise
    '''
    erclist = get_erc(errocodes)
    isin = specific in erclist
    return isin

def case_sensitive_replace(mstring, instring, outstring):
    '''
    replace a sub-string with another sub-string, conserving the case.

    Parameters:
    -----------
    mstring: str
        main string in which replacement to be done
    instring: str
        string that should be replaced with outstring
    outsring: str
        string that replaces the instring

    Returns:
    --------
    string: input string with desired replacements
    '''
    return mstring.replace(instring.upper(), outstring.upper()).replace(
        instring.lower(), outstring.lower())

def datestr2fracday(datestring):
    '''
    transforms data string of form 20180526T030445Z to fractional day of year

    This day of year starts counting on Jan 1st with 1.

    Parameters:
    -----------
    datestring: sring
        date string of the form YYYYMMDDThhmmssZ

    Returns:
    --------
    dof: float
        day of year start count at 1
    '''
    day = datestr2datetime(datestring)
    dof = (day.timetuple().tm_yday + day.timetuple().tm_hour/24. +
           day.timetuple().tm_min/(24.*60))
    return dof
