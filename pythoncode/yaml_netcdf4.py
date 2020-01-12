'''
Module to interchange netCDF4 and yaml formats

This module contains functions to extract the complete info in a netCDF4 file
and save it as a yaml file. Also, a yaml file with the correct structure can be
converted to a nc file.
An nc file can also be simplified and saved as yaml file. This is useful to
create template files.

Author: Martina M. Friedrich

Date: Nov 2019
'''

import copy
import yaml
import netCDF4
import numpy as np
import numpy.ma as ma

try:
    import Tools_MMF as TM
except ImportError:
    from . import Tools_MMF as TM

GASLIST = ["HCHO", "NO2", "H2CO", "SO2", "O3", "BRO"]  # O4 intentionally not
DEFAULT_FILLVALS = {"int32": -2147483647,
                    "float32": np.nan,
                    "float64": np.nan,
                    "uint16": 65535,
                    "int16": -32767,
                    "uint32": 4294967296,
                    "int64": -2147483647
                   }


GEOMS_TYPES = {"REAL": "float32", "DOUBLE": "float64", "STRING": "str"}

def wrap_get_var(mfh, key):
    '''
    needed in order to be able to access also attributes

    Parameters:
    -----------
    mfh: netCDF4 handle
        handle in which to look for the key
    key: str
        key to value to get

    Returns:
    -------
    out: attribute or var
        if var, dict containing info on units etc
    '''
    if "attributes" in key:
        keys = key.split("/attributes/")
        try:
            out = mfh[keys[0]].getncattr(keys[1])
        except AttributeError as exc:
            TM.write_error("ATTR NOT FOUND "+keys[1]+"in"+keys[0], exc,
                           error=False, printtrace=False)
            out = None
    else:
        try:
            out = _get_variables(mfh[key], True, False)
        except KeyError as exc:
            TM.write_error("VAR NOT FOUND "+key+"in"+mfh.path, exc, error=False)
            out = None
        except IndexError as exc:
            TM.write_error("VAR NOT FOUND "+key+"in"+mfh.path, exc, error=False)
            out = None
    return out

def initialize(mfile):
    '''
    Function initialize Tools_MMF with the parameters in the provided file

    Parameters:
    -----------
    mfile: string
        path to configuration file containing general configurations
        frm4doas_tropo.cnf containing at least FILLVALS
    '''
    global DEFAULT_FILLVALS
    _ = TM.initialize(mfile)
    DEFAULT_FILLVALS = {
        "int32": TM.FILLVALS["INT_FILLVAL"],
        "float32": TM.FILLVALS["DBL_FILLVAL"],
        "float64": TM.FILLVALS["DBL_FILLVAL"],
        "uint16": TM.FILLVALS["USHORT_FILLVAL"],
        "int16": TM.FILLVALS["SHORT_FILLVAL"],
        "uint32": TM.FILLVALS["UINT_FILLVAL"],
        "int64": TM.FILLVALS["LONG_FILLVAL"]
        }
    return

def nc_to_dict(netcdf4data, alsodata=False, adjstfv=True):
    '''
    get the structure of a netcdf4 file in a dictionarry structure

    Parameters:
    ----------
    netcdf4data: string
        path to a netCDF4 file from which the structure should be extracted
    alsodata: logical
        if True, also extract the data, if False, extract only the structure
    adjstfv: logical
        if True, adjust FillValues to default ones, if False, do not adjust

    Returns:
    ---------
    mygr: dict
        The dictionarry containing the structure of the netCDF4 file
        all attributs are inside a key "attributes", again as dictionarry
        additionally, all variables have keys for:
        dataype, dimensions, mask, name, units, scale
    '''
    with netCDF4.Dataset(netcdf4data) as myfile:
        mygr = _get_groups(myfile, alsodata, adjstfv)
    return mygr

def yaml_to_nc(yamlfile, ncfile, dimdict=None):
    '''
    turn a yaml file into an nc file. Also need a dict with dimensions

    Parameters:
    -----------
    yamlfile: string
        path to input yaml file
    ncfile: string
        path to the output nc file
    [dimdict: dict
        dictionarry containing as keys all dimensions that are used in the yaml
        file. If not present, yamlfile needs to contain them. ]
    '''
    with open(yamlfile) as mfile:
        if float(yaml.__version__) < 5:
            mdict = yaml.load(mfile)
        else:
            mdict = yaml.full_load(mfile)
    dict_to_nc(mdict, ncfile, dimdict)
    return

def dict_to_nc(mdict, ncfilename, dims=None):
    '''
    turn dictionarry into an nc file. Also need a dict with dimensions

    Parameters:
    -----------
    mdict: dict
        dictionarry containing nc file structure.
    ncfilename: string
        path to the output nc file
    dims: dict
        dictionarry containing as keys all dimensions that are used in the yaml
        file. If dims is None, then the yaml dictionarry for dimensions need to
        be filled

    '''
    #mdict = copy.deepcopy(mdict)
    with netCDF4.Dataset(ncfilename, "w") as ncfile:
        try:
            if dims:
                _make_dimensions(dims, ncfile)
            else:
                _make_dimensions(mdict["dimensions"], ncfile)
                del mdict["dimensions"]
        except KeyError:
            pass  # no dimensions on main level
        for key in mdict.keys():
            if "attributes" in key:
                if "/" in key:  # this means that key is attributes/other
                    ncfile.setncatts({key.split("/")[1]: mdict[key]})
                else:
                    ncfile.setncatts(mdict[key])
            else:
                _create_nc_sub(mdict[key], key, ncfile)
    return

def nc_to_yaml(ncfile, yamlfile, simplyfy=True, alsodata=False, adjustfv=False):
    '''
    extract the structure from a netcdf4 file and write essentials to yaml

    Doublicated groups (doublicated in the sense of same type of group but for
    different gases or different wavelength) and variables are deleted and
    the following replacements (examples) are made:
    123NM --> WAV+NM,  1 --> NUMBER, 123K --> TEMP+K, hcho --> GAS

    Parameters:
    -----------
    ncfile: string
        path to the input netCDF4 file
    yamlfile: string
        path to the output yaml file
    simplyfy: logical
        if True, groups that have the same type of content are not repeated
    alsodata: logical
        if True, also extract the data, if False, extract only structure
    adjustfv: logical
        if True, adjust the FillValues to defaults, if False do not adjust

    The Fillvalues in the ncfile are converted, too. For details, see the
        _get_attributes function
    '''
    mdo = nc_to_dict(ncfile, alsodata, adjustfv)
    if simplyfy:
        md2 = extract_essentials(mdo)
    else:
        md2 = mdo
    with open(yamlfile, "w") as mnf:
        yaml.dump(md2, mnf)
    return

def extract_essentials(mdict_in):
    '''
    some variable names or group names are almost dublicated: detect+delete

    Some groups are replicated and the only difference is e.g. the wavelengt.
    The same is true for some variables. Try to detect these automatically and
    reduce to essentials. Since this is used iteratively, it can also be a
    sub-directory

    Parameters:
    -----------
    mdict_in: dictionarry
        The complete dictionarry represeting the structure of a netCDF file
        This can be either directly from a netCDF file loaded with
        nc_to_dict or also with yaml.load of a saved yaml file created
        with yaml_dump. It can also be a sub-directory.

    Returns:
    --------
    odict: dictionarry
        reduced dictionarry: removed dublicated variables and groups
    '''
    mdict = copy.deepcopy(mdict_in)
    keylist = list(mdict.keys())
    for key in keylist:
        try:
            mdict[key]["path"] = _simplyfy(mdict[key]["path"])
        except (KeyError, TypeError):  # this means it is a group, not a variable
            pass
        except IndexError:
            print("Idx err", key)
        newkey = _simplyfy(key)
        if newkey != key:
            # print (key, newkey)
            mdict[newkey] = mdict[key]
            del mdict[key]
            print("deleted", key, " created ", newkey)
            if key.lower() == key:  # this means it is a variable or attribute
                try:
                    mdict[newkey]["name"] = newkey
                except TypeError:  # this means it is not a variable
                    pass
    for key in mdict:
        if isinstance(mdict[key], dict):
            mdict[key] = extract_essentials(mdict[key])
    return mdict

def _simplyfy(mstr):
    '''
    function to replace a specific number by "str" for e.g QDOAS keys

    If the number was originally followed by "K", then "str" is "TEMP" if it is
    followed by "nm", "str" is "WAV". If followed by nothing, then "NUMBER"

    Parameters:
    -----------
    mstr: str
        string to simplyfy

    Returns:
    --------
    newstr: str
        simplified string
    '''
    if "_" in mstr:
        listkey = mstr.split("_")
        sep = "_"
    elif " " in mstr:
        listkey = mstr.split(" ")
        sep = " "
    else:
        listkey = mstr.split(" ")
        sep = ""
    newstr = []
    for entr in listkey:
        entr = _detect_number(entr)
        entr = _detect_wav(entr)
        entr = _detect_temp(entr)
        entr = _detect_gas(entr)
        newstr.append(entr)
    newstr = sep.join(newstr)
    return newstr

def _get_attributes(parent, adjstfv=False):
    '''
    get all attributes from passed netCDF4 file/ group/ variable into dict

    This does not adjusts the Fillvalue by default. If this should be made,
    the adjustFillValue flag needs to be set to True

    Parameters:
    ----------
    parent: netCDF4 Dataset or group or variable
        netCDF4 entity to get attributes from
    adjstfv: boolean, optional, default: False
        whether or not to adjust the FillValues

    Returns:
    -------
    attr_dict: dictionarry
        dictionarry containing all attributes as key: value pairs
    '''
    attr_dict = {}
    allattr = parent.ncattrs()
    for attr_name in allattr:
        val = parent.getncattr(attr_name)
        if adjstfv and "_FillValue" in attr_name and not isinstance(val, str):
            if isinstance(val, np.float64):
                val = DEFAULT_FILLVALS["float64"]
            elif isinstance(val, np.float32):
                val = DEFAULT_FILLVALS["float32"]
            elif isinstance(val, np.int32):
                val = DEFAULT_FILLVALS["int32"]
            elif isinstance(val, np.int16):
                val = DEFAULT_FILLVALS["int16"]
            elif isinstance(val, np.uint32):
                val = DEFAULT_FILLVALS["uint32"]
            elif isinstance(val, np.uint16):
                val = DEFAULT_FILLVALS["uint16"]
            elif isinstance(val, np.int64):
                val = DEFAULT_FILLVALS["int64"]
            else:
                print(" has unknown dtype ", val, type(val))
        attr_dict[attr_name] = val
    return attr_dict

def _get_groups(parent, alsodata, adjstfv):
    '''
    get all group and variables in passed netCDF4 file or group

    Parameters:
    ----------
    parent: netCDF4 Dataset or group
        netCDF4 entity to get all containing groups and variables from
    alsodata: logical
        should the data also be extracted (True) or just the structure
    adjstfv: lofical
        should fillvalues be adapted (True) to defaults or left as are (False)

    Returns:
    -------
    group_dict: dictionarry
        dictionarry containing all groups/ variables as key: value pairs
        also, the parents attributes are saved under the key "attributes"
    '''
    group_dict = {}
    group_dict["attributes"] = _get_attributes(parent, adjstfv)
    group_dict["dimensions"] = {key: parent.dimensions[key].size for key in
                                parent.dimensions.keys()}
    allgr = parent.groups.keys()
    for gr_name in allgr:
        gr_namen = gr_name.upper()
        group_dict[gr_namen] = _get_groups(parent[gr_name], alsodata, adjstfv)
    allvars = parent.variables.keys()
    for var_name in allvars:
        group_dict[var_name] = _get_variables(
            parent[var_name], alsodata, adjstfv)
    return group_dict

def _get_variables(parent, alsodata, adjstfv):
    '''
    get a variable in passed netCDF4 file or group

    Parameters:
    ----------
    parent: netCDF4 variable
        netCDF4 entity to get all containing variables from
    alsodata: logical
        should the data also be extracted (True) or just the structure
    adjstfv: lofical
        should fillvalues be adapted (True) to defaults or left as are (False)

    Returns:
    -------
    var_dict: dictionarry
        dictionarry containing all info about a variable, if alsodata=True
        also contains the data
    '''
    var_dict = {}
    var_dict["attributes"] = _get_attributes(parent, adjstfv)
    dtp = str(parent.datatype)
    if "string" in dtp:
        dtp = "str"
    var_dict["datatype"] = str(dtp)
    var_dict["dimensions"] = (tuple(parent.dimensions))
    var_dict["mask"] = str(parent.mask)
    var_dict["name"] = parent.name
    var_dict["path"] = parent.group().path
    var_dict["scale"] = str(parent.scale)
    if alsodata:
        var_dict["value"] = parent[:].astype(dtp)
        if adjstfv:
            newfillval = var_dict["attributes"]["_FillValue"]
            oldfillval = parent._FillValue
            if "str" not in var_dict["datatype"]:
                if np.isnan(oldfillval):
                    try:
                        var_dict["value"][np.isnan(var_dict["value"].data)] = newfillval
                    except:
                        var_dict["value"][np.isnan(var_dict["value"])] = newfillval
                else:
                    try:
                        var_dict["value"][var_dict["value"].data == oldfillval] = newfillval
                    except:
                        var_dict["value"][var_dict["value"] == oldfillval] = newfillval
                if isinstance(var_dict["value"], ma.core.MaskedArray):
                    var_dict["value"].fill_value = newfillval
                    if isinstance(var_dict["value"].mask, np.bool_):
                        var_dict["value"].mask = False
                    if np.isnan(newfillval):
                        var_dict["value"].mask[np.isnan(var_dict["value"])] = True
                    else:
                        var_dict["value"].mask[var_dict["value"] == newfillval] = True
                else:
                    var_dict["value"] = ma.array(
                        var_dict["value"], fill_value=newfillval)
    else:
        var_dict["value"] = ""
    return var_dict

def _detect_number(entr):
    '''
    replace plain integer by the word NUMBER

    Parameters:
    -----------
    entr: string
        string to be possibly replaced with the word NUMBER

    Returns:
    --------
    entr: string
        either the original string, or, if the string could be converted to int
        the word NUMBER is returned
    '''
    try:
        entr = int(entr)
        entr = "NUMBER"
    except ValueError:
        pass
    return entr

def _detect_wav(mstr):
    '''
    replace a number followed by "NM" by the string WAV+NM

    Parameters:
    -----------
    entr: string
        string to be possibly replaced with WAV+NM

    Returns:
    --------
    val: string
        either the original string, or, if the string could be converted to
        a float after stripping "NM" on the right, WAV+NM
    '''
    try:
        _ = float(mstr.upper().rstrip("NM"))
        val = "WAV+NM"
    except ValueError:
        val = mstr
    return val

def _detect_temp(mstr):
    '''
    replace a number followed by "K" by the string TEMP+K

    Parameters:
    -----------
    entr: string
        string to be possibly replaced with TEMP+K

    Returns:
    --------
    entr: string
        either the original string, or, if the string could be converted to
        a float after stripping "K" on the right, TEMP+K
    '''
    try:
        _ = float(mstr.upper().rstrip("K"))
        val = "TEMP+K"
    except ValueError:
        val = mstr
    return val

def _detect_gas(mstr):
    '''
    check whether a string is in GASLIST. If so, return GAS, otherwise string

    Parameters:
    -----------
    mstr: string
        string to test

    Returns:
    --------
    string
        either the original string or GAS
    '''
    if mstr.upper() in GASLIST:
        return "GAS"
    else:
        return mstr

def _make_dimensions(dims, parent):
    '''
    create dimensions in nc file or group

    Parameters:
    -----------
    dims: dict
        a dictionarry of dimensions that need to be created in parent
    parent: nc group or file handle
        where to create the dimensions
    '''
    for dim in dims.keys():
        parent.createDimension(dim, int(dims[dim]))
    return

def _create_nc_sub(mdict, key, parent):
    '''
    create group or variable (from mdict) inside nc handle parent under key

    Parameters:
    -----------
    mdict: dict
        dictionarry what to create
    key: string
        name of the group/ variable to create
    parent: nc file/ group handle
        where to create this group or variable
    '''
    if mdict is None:
        print("None?", key)
        return
    #if key.upper() == key:
    if key == "attributes":
        return
    if isinstance(mdict, dict) and "value" not in mdict.keys():
        child = parent.createGroup(key)
        try:
            child.setncatts(mdict["attributes"])
        except (KeyError, TypeError):
            noattr = True
            for skey in mdict.keys():
                if "attributes" in skey:
                    noattr = False
                    try:
                        attrkey = skey.split("/")[1]
                        child.setncatts({attrkey: mdict[skey]})
                    except IndexError:
                        print("not dict, no entries: ", skey)
            if noattr:
                pass
                #print("no attributes for ", key)
            #attrd = {}
            #for entr in ["description", "units", "_FillValue"]:
            #    try:
            #        attrd[entr] = mdict[entr]
            #    except KeyError:
            #        print("failed to add some attributes:")
            #        print("  ", key, entr)
        try:
            _make_dimensions(mdict["dimensions"], child)
            del mdict["dimensions"]
        except:
            pass
        try:
            if isinstance(mdict, dict):
                for key2 in mdict.keys():
                    _create_nc_sub(mdict[key2], key2, child)
            else:
                print("is type: ", type(mdict), " for ", key)
        except AttributeError as atr:
            print("problem: ", atr)
            print("not a group with content: ", key)
            print("caused by : ", key2)
            print("type: ", type(mdict))
    else:
        if isinstance(mdict, dict) and "attributes" not in key:
            try:
                attrd = copy.deepcopy(mdict["attributes"])
            except KeyError:
                attrd = {}
                for skey in mdict.keys():
                    if "attribute" in skey:
                        name = skey.split("/")[1]
                        attrd[name] = mdict[skey]
            for entr in ["description", "units", "_FillValue", "source"]:
                try:
                    attrd[entr] = mdict[entr]
                except KeyError:
                    pass
            if "_FillValue" in attrd.keys():
                fillval = attrd["_FillValue"]
            else:
                fillval = DEFAULT_FILLVALS[mdict["datatype"]]
            if not isinstance(mdict["dimensions"], str):
                dimensions = tuple(d for d in mdict["dimensions"])
            else:
                dimensions = (mdict["dimensions"],)
            try:
                newvar = parent.createVariable(
                    mdict["name"], mdict["datatype"],
                    dimensions, fill_value=fillval)
            except RuntimeError as err:
                print(err)
            except KeyError:  # this is for geoms
                newvar = parent.createVariable(
                    mdict["name"],
                    GEOMS_TYPES[mdict["attributes"]["VAR_DATA_TYPE"]],
                    dimensions, fill_value=mdict["_FillValue"])
            try:
                del attrd["_FillValue"]
            except KeyError:
                print("no fillval for", key)
            try:
                newvar.setncatts(attrd)
            except Exception as exc:
                print("problems: ", mdict["name"], " attr problem: ", exc)
                for atr in attrd:
                    print("    ", atr, attrd[atr])
            try:
                #if mdict["name"].lower() in "air_mass_factor_of_o3":
                #    import pdb
                #    pdb.set_trace()
                newvar[:] = mdict["value"]
            except IndexError:
                print("problem: ", mdict["name"], len(mdict["value"]), mdict["dimensions"])
                newvar[:] = np.array([mdict["value"]])
            except ValueError as err:
                print(" has no value", key)
                print("err:", err)
            except TypeError as err:
                print("has a problem", key)
                print(err)
                TM.write_error("ee", err)
                print(mdict.keys(), mdict["value"])
                #newvar[:] = [mdict["value"]]

    return

def collapse_keys(mobj, last_key=None, mdc=None):
    '''
    nested dict to dict with 1 level and concatenated keys

    translate a nested list of dictionarries into one with keys joind by "/"
    the variables are still dictionarries. However, there might be attributes
    in the groups. This might cause problems and need extra checks if it is
    assumed that the last key not followed by any / is a variable.

    Parameters:
    -----------
    mobj: dict
        dict that has nested dicts
    [last_key: str
        if this is already a nested dict, this is the key of that dict]
    [mdc: dict
        if this is already nested dict, this is the previous dict]

    Returns:
    -------
    mdc: dict
        flat dictionarry with collapsed keys containing "/" for each level

    Example:
    --------
    mdict["one"]["one_1"] = "a"
    mdict["one"]["one_2"] = "b"
    mdict["two"]["two_1"] = "c"
    -->
    newdict["one/one_1"] = "a"
    newdict["one/one_2"] = "b"
    newdict["two/two_1"] = "c"
    '''
    if mdc is None:
        mdc = {}
    if isinstance(mobj, dict) and len(mobj.keys()) > 0 and "value" not in mobj.keys():
        for key in mobj.keys():
            if last_key is not None:
                last_key_h = "/".join([last_key, key])
            else:
                last_key_h = key
            #if key not in ["dimensions", "attributes"]:
            md2 = collapse_keys(mobj[key], last_key_h, mdc)
            md2.update(md2)
        mdc.update(md2)
    else:
        mdc[last_key] = mobj
    return mdc

def myupdatedict(dict1, dict2):
    '''
    make a recurse dictionarry update: update dict1 with dict2

    Parameters:
    ----------
    dict1: dict
        dictionarry to update
    dict2: dict
        dictionarry used for the update
    '''
    newdict = {}
    for key in dict1:
        if key not in dict2.keys():
            newdict[key] = dict1[key]
        else:
            if isinstance(dict1[key], dict) and isinstance(dict2, dict):
                newdict[key] = myupdatedict(dict1[key], dict2[key])
            else: # if they are not both dict, take dict2
                newdict[key] = dict2[key]
    for key in dict2:
        if key not in dict1.keys():
            newdict[key] = dict2[key]
    return newdict

def expand_keys(mdict):
    '''
    make nested dictionarries of 1 layer dict where keys are joint with /

    This is the opposite of collapse_keys

    Parameters:
    ---------
    mdict: dict
        one-layer dictionarry containing the same info as mdict

    Returns:
    --------
    newdict: dict
        nested dictionarry

    '''
    newdict = {}
    doagain = False
    for key in mdict.keys():
        key = key.lstrip("/")
        keylist = key.split("/")
        lastkey = keylist[-1]
        restkey = "/".join(keylist[:-1])
        if restkey:
            if restkey not in newdict.keys():
                newdict[restkey] = {lastkey: mdict[key]}
            else:
                if lastkey in newdict[restkey].keys():
                    newdict[restkey][lastkey] = myupdatedict(
                        newdict[restkey][lastkey], mdict[key])
                else:
                    newdict[restkey][lastkey] = mdict[key]
        else:
            if key not in newdict.keys():
                newdict[key] = mdict[key]
            else:
                newdict[key] = myupdatedict(newdict[key], mdict[key])
    for key in newdict:
        if "/" in key:
            doagain = True
    if doagain:
        newdict = expand_keys(newdict)
    return newdict

# this is just for testing stuff
# class ndarray(np.ma.core.MaskedArray):
#     def __new__(cls, mdict):
#        obj = np.ma.core.MaskedArray(mdict).view(cls)
#        obj.extra_info = "nop"
#        return obj
#
# def ndarray(loader, node):
#    '''
#    enables ndarray types in yaml
#    '''
#    mapping = loader.construct_sequence(node, deep=True)
#    mat = np.array(mapping)
#    return mat
#
# yaml.add_constructor(u"tag:yaml.org,2002:ndarray", ndarray)
