"""Set of common functions that are used in loads of scripts."""


import re
import ROOT
import os
from subprocess import call
import sys
from sys import platform as _platform
import numpy as np
import math
import argparse
from array import array
from itertools import chain
from collections import OrderedDict
import pickle

np.set_printoptions(edgeitems=3, infstr='Infinity',
                    linewidth=220, nanstr='nan', precision=6,
                    suppress=False, threshold=sys.maxsize, formatter=None)


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.TH1.SetDefaultSumw2(True)


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    """Make argparse respect space formatting (newlines, etc) AND show defaults"""
    pass


def open_pdf(pdf_filename):
    """Open a PDF file using system's default PDF viewer."""
    if _platform.startswith("linux"):
        # linux
        call(["xdg-open", pdf_filename])
    elif _platform == "darwin":
        # OS X
        call(["open", pdf_filename])
    elif _platform == "win32":
        # Windows
        call(["start", pdf_filename])


#
# Filepath/directory fns
#
def cleanup_filepath(filepath):
    """Resolve any sym links, env vars, ~, etc, and return absolute path."""
    thing = os.path.abspath(os.path.expandvars(os.path.expanduser(filepath)))
    if os.path.islink(thing):
        # because otherwise readlink throws if not a link
        linkpath = os.readlink(thing)
        if not os.path.isabs(linkpath):
            # convert to abs path, otherwise it uses the relpath from the py
            # file location (which will be wrong)
            dirname = os.path.dirname(thing)
            newpath = os.path.join(dirname, linkpath)
            return newpath
        else:
            return linkpath
    else:
        return thing


def get_full_path(filepath):
    """Return absolute directory of filepath.
    Resolve any environment vars, ~, sym links(?)"""
    return os.path.dirname(cleanup_filepath(filepath))


def check_file_exists(filepath):
    """Check if file exists. Can do absolute or relative file paths."""
    return os.path.isfile(cleanup_filepath(filepath))


def check_dir_exists(filepath):
    """Check if directory exists."""
    return os.path.isdir(cleanup_filepath(filepath))


def check_dir_exists_create(filepath):
    """Check if directory exists. If not, create it."""
    if not check_dir_exists(filepath):
        os.makedirs(cleanup_filepath(filepath))


#
# ROOT specific fns, like opening files safely
#


class TFileCacher(object):
    """Simple class to cache TFile & objects already "gotten" from TFile

    Can either do mytfile.get(x) or just mytfile[x]

    TODO: Should this just inherit from TFile instead?
    """
    def __init__(self, filename):
        self.tfile = open_root_file(filename)
        self.objects = {}

    def get(self, obj_name):
        return self[obj_name]

    def Get(self, obj_name):
        # uppercase as well, like in TFile
        return self[obj_name]

    def __getitem__(self, obj_name):
        # really hoping all these have unique names
        if obj_name not in self.objects:
            self.objects[obj_name] = get_from_tfile(self.tfile, obj_name)
        return self.objects[obj_name]

    def Close(self):
        self.tfile.Close()


def open_root_file(filename, mode="READ"):
    """Safe way to open ROOT file. Could be improved."""
    clean_filename = cleanup_filepath(filename)
    if mode in ["READ", "UPDATE"]:
        if not check_file_exists(clean_filename):
            raise IOError("No such file %s" % clean_filename)
    f = ROOT.TFile(clean_filename, mode)
    if f.IsZombie() or not f:
        raise IOError("Can't open TFile %s" % clean_filename)
    return f


def exists_in_file(tfile, obj_name):
    """Check if object exists in TFile.

    Also handles directory structure, e.g. histograms/central/pt_1
    """
    parts = obj_name.split("/")
    current_obj = tfile
    for p in parts:
        if current_obj.GetListOfKeys().Contains(p):
            current_obj = current_obj.Get(p)
        else:
            return False
    return True


def get_from_tfile(tfile, obj_name, info=False):
    """Get some object from ROOT TFile with checks."""
    if info:
        print("Getting %s" % obj_name)
    obj = tfile.Get(obj_name)
    if obj == None:
        raise IOError("No object named %s in %s" % (obj_name, tfile.GetName()))
    else:
        if isinstance(obj, (ROOT.TH1, ROOT.TH2, ROOT.TH3)):
            obj.SetDirectory(0)  # Free it from the TFile ownership
        return obj


def grab_obj_from_file(file_name, obj_name):
    """Get object names obj_name from ROOT file file_name"""
    input_file = open_root_file(file_name)
    obj = get_from_tfile(input_file, obj_name)
    if isinstance(obj, (ROOT.TH1, ROOT.TH2)):
        # THIS ORDER IS VERY IMPORTANT TO AVOID MEMORY LEAKS
        new_obj = obj.Clone(get_unique_str())
        new_obj.SetDirectory(0)  # Ownership kludge
        input_file.Close()
        return new_obj
    else:
        return obj


def check_root_obj(obj):
    """Check ROOT object is fine"""
    if obj.IsZombie() or obj is None:
        raise IOError("object is unavailable")


def close_tfile(obj):
    """Close TFile, after checking it's open"""
    this_obj = obj
    if isinstance(obj, TFileCacher):
        this_obj = obj.tfile
    if isinstance(this_obj, ROOT.TFile) and this_obj.IsOpen():
        this_obj.Close()


def check_exp(n):
    """
    Checks if number has stupidly larger exponent

    Can occur is using buffers - it just fills unused bins with crap
    """

    from math import fabs, log10, frexp
    m, e = frexp(n)
    return fabs(log10(pow(2, e))) < 10


def get_xy(graph):
    """
    Return lists of x, y points from a graph, because it's such a PITA
    """
    xarr = list(np.ndarray(graph.GetN(), 'd', graph.GetX()))
    yarr = list(np.ndarray(graph.GetN(), 'd', graph.GetY()))
    return xarr, yarr


def get_exey(graph):
    """
    Return lists of errors on x, y points from a graph, because it's such a PITA
    """
    xarr = list(np.ndarray(graph.GetN(), 'd', graph.GetEX()))
    yarr = list(np.ndarray(graph.GetN(), 'd', graph.GetEY()))
    return xarr, yarr


def th2_to_arr(h, do_errors=False):
    """Convert TH2 to 2D numpy array"""
    arr = np.zeros((h.GetNbinsY(), h.GetNbinsX()), dtype=float)
    err_arr = np.zeros((h.GetNbinsY(), h.GetNbinsX()), dtype=float)
    for x_ind in range(1, h.GetNbinsX() + 1):
        for y_ind in range(1, h.GetNbinsY() + 1):
            # arr[x_ind-1][y_ind-1] = h.GetBinContent(x_ind, y_ind)
            arr[y_ind-1][x_ind-1] = h.GetBinContent(x_ind, y_ind)
            if do_errors:
                # err_arr[x_ind-1][y_ind-1] = h.GetBinError(x_ind, y_ind)
                err_arr[y_ind-1][x_ind-1] = h.GetBinError(x_ind, y_ind)
    return arr, err_arr


def make_normalised_TH2(hist, norm_axis, recolour=True, do_errors=False):
    norm_axis = norm_axis.upper()
    possible_opts = ['X', 'Y']
    if norm_axis not in possible_opts:
        raise RuntimeError("norm_axis must be one of %s" % possible_opts)
    norm_axis = norm_axis.upper()

    # easiest way to cope with x or y is to just get a 2D matrix of values,
    # can then do transpose if necessary
    arr, err_arr = th2_to_arr(hist, do_errors)

    if norm_axis == 'X':
        arr = arr.T
        err_arr = err_arr.T

     # all operations here are done per row (i.e. per y bin)
    if recolour:
        # can set so the maximum in each bin is the same,
        # scale other bins accordingly
        # this retain the colour scheme for each set of bins
        maxes = arr.max(axis=1)  # maximum per row
        maxesT = maxes[:, None]
        arr = np.divide(arr, maxesT, where=maxesT!=0, out=np.zeros_like(arr))
        err_arr = np.divide(err_arr, maxesT, where=maxesT!=0, out=np.zeros_like(arr))
        # for ind, xbin in enumerate(arr):
        #     if xbin.max() > 0:
        #         factor = xbin.max()
        #         arr[ind] = xbin / factor
        #         err_arr[ind] = err_arr[ind] / factor
    else:
        # alternatively, can rescale so sum over bins = 1
        for ind, xbin in enumerate(arr):
            summed = arr.sum(axis=1)
            summedT = summed[:, None]
            arr = np.divide(arr, summedT, where=summedT!=0, out=np.zeros_like(arr))
            err_arr = np.divide(err_arr, summedT, where=summedT!=0, out=np.zeros_like(arr))
            # if xbin.sum() != 0:
            #     factor = xbin.sum()
            #     arr[ind] = xbin / factor
            #     err_arr[ind] /= factor

    if norm_axis == 'X':
        arr = arr.T
        err_arr = err_arr.T

    # Create new hist object - MUST do it this way to get Z range correct
    new_histname = hist.GetName() + "_norm" + norm_axis
    # hnew = ROOT.TH2F(hist)  # TODO: determine if TH2F, TH2D...
    if type(hist) == ROOT.TH2F:
        hnew = ROOT.TH2F(hist)  # TODO: determine if TH2F, TH2D...
    elif type(hist) == ROOT.TH2D:
        hnew = ROOT.TH2D(hist)  # TODO: determine if TH2F, TH2D...
    else:
        raise RuntimeError("Unknown 2D hist type")
    hnew.SetName(new_histname)
    for y_ind, y_arr in enumerate(arr, 1):
        for x_ind, val in enumerate(y_arr, 1):
            hnew.SetBinContent(x_ind, y_ind, val)
            if do_errors:
                hnew.SetBinError(x_ind, y_ind, err_arr[y_ind-1][x_ind-1])
#     hnew.SetAxisRange(0.5, 1., 'Z')
    return hnew


def th1_to_arr(hist):
    return np.array([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX()+1)]), np.array([hist.GetBinError(i) for i in range(1, hist.GetNbinsX()+1)])


def get_list_of_element_names(thing):
    """Get list of key names for given thing in a ROOT file"""
    return [x.GetName() for x in thing.GetListOfKeys()]


def get_list_of_objects_in_dir(filename, dirname):
    """Get list of elements in a TDirectory"""
    f = open_root_file(filename)  # need to keep this in memory otherwise the get_list... segfaults
    d = get_from_tfile(f, dirname)
    return get_list_of_element_names(d)


def get_unique_str():
    return ROOT.TUUID().AsString()


def thstack_to_th1(hst):
    """Add all hists in a THStack into one TH1"""
    htotal = None
    hists = hst.GetHists()
    for h in hists:
        if htotal is None:
            htotal = h.Clone(get_unique_str())
        else:
            htotal.Add(h)
    return htotal


def get_hist_mean_rel_error(hist):
    """Get average relative error from histogram bins (i.e. error / contents)

    Parameters
    ----------
    hist : ROOT.TH1 (or descendents)
        Description

    Returns
    -------
    float
        average relative error
    """
    nbins = hist.GetNbinsX()
    contents = np.array([hist.GetBinContent(i) for i in range(1, nbins+1)])
    errors = np.array([hist.GetBinError(i) for i in range(1, nbins+1)])
    # only care about bins with something in them
    mask = contents > 0
    # if no bins have any contents, just return 0
    if not mask.any():
        return 0.
    rel_err = np.divide(errors[mask], contents[mask])
    rel_err[np.isinf(rel_err)] = 0.
    rel_err[np.isnan(rel_err)] = 0.
    return np.average(rel_err)


def get_jet_config_from_dirname(dirname):
    """Get jet configuration from string

    Parameters
    ----------
    dirname : str
        Description

    Returns
    -------
    (str, str)
        jet radius and PUS

    Raises
    ------
    RuntimeError
        If it can't find it
    """

    m = re.search(r'ak([\d]+)(chs|puppi)', dirname.lower())
    if m and len(m.groups()) == 2:
        return m.groups()
    else:
        print("dirname:", dirname)
        if m:
            print(m)
        raise RuntimeError("cannot determine jet config from dirname %s" % dirname)


def get_bin_edges(hist, axis):
    """Get array of bin edges from hist. Must specify which axis to use."""
    axis = axis.lower()
    if axis not in ['x', 'y']:
        raise RuntimeError("get_bin_edges axis must be x or y")
    ax = hist.GetXaxis() if axis == "x" else hist.GetYaxis()
    bins = [ax.GetBinLowEdge(i) for i in range(1, ax.GetNbins()+2)]
    return bins


class Marker(object):

    shape_dict = OrderedDict()
    shape_dict['circle'] = {'filled': 20, 'open': 24}
    shape_dict['square'] = {'filled': 21, 'open': 25}
    shape_dict['triangleUp'] = {'filled': 22, 'open': 26}
    shape_dict['triangleDown'] = {'filled': 23, 'open': 32}
    shape_dict['star'] = {'filled': 29, 'open': 30}
    shape_dict['diamond'] = {'filled': 33, 'open': 27}
    shape_dict['cross'] = {'filled': 34, 'open': 28}
    shape_dict['crossX'] = {'filled': 47, 'open': 46}
    shape_dict['doubleDiamond'] = {'filled': 43, 'open': 42}
    shape_dict['3star'] = {'filled': 39, 'open': 37}


    def __init__(self, shape='circle', filled=True):
        self.fill_state = filled
        self.shape_state = shape

    @staticmethod
    def get(shape, filled=True):
        if shape not in Marker.shape_dict.keys():
            raise RuntimeError("Unknown marker shape %s" % (shape))
        return Marker.shape_dict[shape]['filled' if filled else 'open']

    def cycle(self, filled=True, cycle_filling=False, only_cycle_filling=False):
        if only_cycle_filling:
            self.fill_state = filled
            while True:
                yield self.get(self.shape_state, self.fill_state)
                self.fill_state = not self.fill_state
        else:
            # start cycle at self.shape_state
            keys = list(Marker.shape_dict.keys())
            ind = keys.index(self.shape_state)
            keys = keys[ind:] + keys[:ind]
            for k in chain(keys):
                self.fill_state = filled
                yield self.get(k, self.fill_state)
                if cycle_filling:
                    self.fill_state = not self.fill_state
                    yield self.get(k, self.fill_state)


def set_french_flag_palette():
    """Custom colour scheme - french flag colouring"""
    NRGBs = 5
    ff_ncont = 256
    # Set max & min such that 0 is in the middle
    stops = [ 0.00, 0.49, 0.51, 0.52, 1.00 ]
    red = [ 0.00, 1.00, 1.00, 1.00, 1.00]
    green = [ 0.00, 1.00, 1.00, 1.00, 0.00 ]
    blue = [ 1.00, 1.00, 1.00, 1.00, 0.00 ]
    stopsArray = array('d', stops)
    redArray = array('d', red)
    greenArray = array('d', green)
    blueArray = array('d', blue)
    french_flag_cols = ROOT.TColor.CreateGradientColorTable(NRGBs, stopsArray, redArray, greenArray, blueArray, ff_ncont)
    french_flag_palette = array('i', [french_flag_cols+i for i in range(ff_ncont)])
    ROOT.gStyle.SetPalette(ff_ncont, french_flag_palette)


def set_log_french_flag_palette():
    """Custom colour scheme - french flag colouring but colour small values non-white"""
    NRGBs = 5
    ff_ncont = 512
    # Set max & min such that 0 is in the middle
    delta = 1E-20
    # TODO: actually have log colour scale?
    stops = [ 0.00, 0.5-(delta), 0.5, 0.5+(delta), 1.00 ]
    # red =   [ 0.00, 0.00, 0.50, 0.98, 1.00]
    # green = [ 0.00, 1.00, 0.93, 0.86, 0.00 ]
    # blue =  [ 1.00, 1.00, 0.98, 0.96, 0.00 ]
    red =   [ 0.00, 0.00, 1.00, 0.98, 1.00]
    green = [ 0.00, 1.00, 1.00, 0.86, 0.00 ]
    blue =  [ 1.00, 1.00, 1.00, 0.96, 0.00 ]
    stopsArray = array('d', stops)
    redArray = array('d', red)
    greenArray = array('d', green)
    blueArray = array('d', blue)
    french_flag_cols = ROOT.TColor.CreateGradientColorTable(NRGBs, stopsArray, redArray, greenArray, blueArray, ff_ncont)
    french_flag_palette = array('i', [french_flag_cols+i for i in range(ff_ncont)])
    ROOT.gStyle.SetPalette(ff_ncont, french_flag_palette)


def symmetrize_h2d_z_limits(h2d):
    """Set TH2 z axis to have same + & - limits - good for use with French flag palette"""
    max_z = max(h2d.GetMaximum(), abs(h2d.GetMinimum()))
    h2d.SetMaximum(max_z)
    h2d.SetMinimum(-max_z)


def get_colour_seq(index, total_length):
    """Get colour from current palette that corresponds to (index / total_length) fraction in the palette"""
    num_colours = ROOT.TColor.GetPalette().fN
    return ROOT.TColor.GetColorPalette(int(num_colours * index / total_length))


# Various helper functions
# ------------------------------------------------------------------------------

SPACE_REPLACEMENT_CHAR = "-"


def no_space_str(s):
    return s.replace(" ", SPACE_REPLACEMENT_CHAR) if s else s


def str_restore_space(s):
    return s.replace(SPACE_REPLACEMENT_CHAR, " ") if s else s


def get_dict_item_sizes(this_dict, recursive=True):
    """Get size of elements in dict.

    Can do dicts recursively

    Uses pickle to get the "true" size
    """
    size_dict = {}
    # print("size dict:", this_dict)
    for k, v in this_dict.items():
        # print("Getting...", k, type(v), v)
        this_key = "[%s]" % k
        if isinstance(v, dict) and recursive:
            tsd = get_dict_item_sizes(v, recursive)
            # print("k,tsd:", k, tsd)
            size_dict.update({"%s%s" % (this_key, kk) : vv for kk, vv in tsd.items()})
            # print("new size_dict from dict", size_dict)
        else:
            if v is None:
                continue
            # print(k, v, type(v))
            obj_pkl = pickle.dumps(v)
            size_dict[this_key] = sys.getsizeof(obj_pkl)
            # print("new size_dict", size_dict)

    return size_dict


def print_dict_item_sizes(this_dict, descending=True, recursive=True):
    """Print out size of elements in dict

    If descending is True, then do largest size first.

    If recursive is True, then do dicts recursively.
    """
    size_dict = get_dict_item_sizes(this_dict, recursive)
    # FIXME broken in python2?
    # print out, sorted by size
    sorted_dict = {k:v for k,v in sorted(size_dict.items(), key=lambda x: x[1], reverse=descending)}
    for k, v in sorted_dict.items():
        print(k, v)


def nsf(num, n=1):
    """Get n-Significant figures from float"""
    numstr = ("{0:.%ie}" % (n-1)).format(num)
    return float(numstr)


# Various methods to convert between ROOT things and numpy
# ------------------------------------------------------------------------------

def tmatrixdsparse_to_ndarray(matrix):
    ndarr = np.zeros(shape=(matrix.GetNrows(), matrix.GetNcols()))

    rows_A = matrix.GetRowIndexArray()
    cols_A = matrix.GetColIndexArray()
    data_A = matrix.GetMatrixArray()
    for iy in range(matrix.GetNrows()):
        for indexA in range(rows_A[iy], rows_A[iy+1]):
            ix = cols_A[indexA]
            # print([x for x in self.GetXToHist()])
            # TODO: care about orientation?
            ndarr[iy, ix] = data_A[indexA]
    return ndarr


def th2_to_tmatrixd(hist, include_uflow=False, include_oflow=False):
    n_rows = hist.GetNbinsY()
    n_cols = hist.GetNbinsX()

    # ignore for now as too complicated
    # if include_uflow:
    #     n_rows += 1
    #     n_cols += 1
    # if include_oflow:
    #     n_rows += 1
    #     n_cols += 1

    # taken from https://root.cern.ch/doc/master/TH2_8cxx_source.html#l03739
    m = ROOT.TMatrixD(n_rows, n_cols)
    ilow = m.GetRowLwb()
    iup  = m.GetRowUpb()
    jlow = m.GetColLwb()
    jup  = m.GetColUpb()
    for i in range(ilow, iup+1):
        for j in range(jlow, jup+1):
            m[i,j] = hist.GetBinContent(j-jlow+1,i-ilow+1)
    return m


def get_th1_bin_centers(h):
    # TODO maintain same shape(1, n) as in th1_to_array?
    centers = np.array([h.GetBinLowEdge(i) + 0.5*h.GetBinWidth(i) for i in range(1, h.GetNbinsX()+1)])
    return centers


def get_th1_bin_widths(h):
    # TODO maintain same shape(1, n) as in th1_to_array?
    widths = np.array([h.GetBinWidth(i) for i in range(1, h.GetNbinsX()+1)])
    return widths


def th1_to_ndarray(hist_A, oflow_x=False):
    """Convert TH1 to numpy ndarray"""
    ncol = hist_A.GetNbinsX()
    if oflow_x:
        ncol += 2
    contents = np.zeros(shape=(1, ncol), dtype=np.float64)
    errors = np.zeros(shape=(1, ncol), dtype=np.float64)

    # Get ROOT indices to loop over
    x_start = 0 if oflow_x else 1
    x_end = hist_A.GetNbinsX()
    if oflow_x:
        x_end += 1

    # x_ind for numpy as always starts at 0
    # ix for ROOT
    for x_ind, ix in enumerate(range(x_start, x_end+1)):
        contents[0][x_ind] = hist_A.GetBinContent(ix)
        errors[0][x_ind] = hist_A.GetBinError(ix)

    # check sparsity
    return contents, errors


def ndarray_to_th1(nd_array, has_oflow_x=False, offset=0.5, bins=None):
    """Convert numpy ndarray row vector to TH1, with shape (1, nbins)

    Use has_oflow_x to include the under/overflow bins
    """
    nbinsx = nd_array.shape[1]
    nbins_hist = nbinsx
    if has_oflow_x:
        nbins_hist -= 2

    # need the 0.5 offset to match TUnfold
    if bins == None:
        h = ROOT.TH1F(get_unique_str(), "", nbins_hist, offset, nbins_hist+offset)
    else:
        if len(bins) != nbinsx+1:
            raise IndexError("len(bins) != nbinsx + 1")
        h = ROOT.TH1F(get_unique_str(), "", nbins_hist, bins)

    x_start = 1
    x_end = nbins_hist

    if has_oflow_x:
        x_start = 0
        x_end = nbins_hist+1

    for x_ind, ix in enumerate(range(x_start, x_end+1)):
        h.SetBinContent(ix, nd_array[0][x_ind])
        h.SetBinError(ix, math.sqrt(abs(nd_array[0][x_ind])))
        #FIXME how to do errors
    return h


def th2_to_ndarray(hist_A, oflow_x=False, oflow_y=False):
    """Convert TH2 to numpy ndarray

    Don't use verison in common_utils - wrong axes?
    """
    ncol = hist_A.GetNbinsX()
    if oflow_x:
        ncol += 2
    nrow = hist_A.GetNbinsY()
    if oflow_y:
        nrow += 2

    contents = np.zeros(shape=(nrow, ncol), dtype=np.float64)
    errors = np.zeros(shape=(nrow, ncol), dtype=np.float64)
    # access via contents[irow][icol]

    # Get ROOT indices to loop over
    y_start = 0 if oflow_y else 1
    y_end = hist_A.GetNbinsY()
    if oflow_y:
        y_end += 1

    x_start = 0 if oflow_x else 1
    x_end = hist_A.GetNbinsX()
    if oflow_x:
        x_end += 1

    # y_ind, x_ind for numpy as always starts at 0
    # iy, ix for ROOT
    for y_ind, iy in enumerate(range(y_start, y_end+1)):
        for x_ind, ix in enumerate(range(x_start, x_end+1)):
            contents[y_ind][x_ind] = hist_A.GetBinContent(ix, iy)
            errors[y_ind][x_ind] = hist_A.GetBinError(ix, iy)

    return contents, errors


def ndarray_to_th2(data, offset=0, binsx=None, binsy=None):
    nbinsy, nbinsx = data.shape
    bins_x = array('d', [x+offset for x in range(1, nbinsx+2)]) # e.g. offset = -0.5 for TUnfold
    bins_y = array('d', [x+offset for x in range(1, nbinsy+2)])
    if binsx is not None:
        if len(binsx) != nbinsx+1:
            raise IndexError("binsx wrong size")
        bins_x = binsx
    if binsy is not None:
        if len(binsy) != nbinsy+1:
            raise IndexError("binsy wrong size")
        bins_y = binsy
    h = ROOT.TH2D(get_unique_str(), "", nbinsx, bins_x, nbinsy, bins_y)
    for ix in range(nbinsx):
        for iy in range(nbinsy):
            h.SetBinContent(ix+1, iy+1, data[iy,ix])
            h.SetBinError(ix+1, iy+1, 0)
    return h


def normalise_ndarray(matrix, by):
    if by == 'col':
        matrix = matrix.T # makes life a bit easier
    for i in range(matrix.shape[0]):
        row_sum = matrix[i].sum()
        if row_sum != 0:
            matrix[i] = matrix[i] / row_sum
    if by == 'col':
        return matrix.T
    else:
        return matrix


def shift_to_covariance(hist):
    """Convert 1D shift histogram to 2D covariance matrix, using V = x x^T"""
    xax = hist.GetXaxis()
    bins = array('d', [xax.GetBinLowEdge(i) for i in range(1, hist.GetNbinsX()+2)])
    nbins = len(bins) - 1
    h2d = ROOT.TH2D("covariance_" + hist.GetName(), "Covariance;%s;%s" % (xax.GetTitle(), xax.GetTitle()), nbins, 0, nbins, nbins, 0, nbins)
    values = np.array([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX()+1)])
    values = values.reshape(len(values), 1)  # turn into column vector
    cov_values = values.dot(values.T)
    for ix in range(nbins):
        for iy in range(nbins):
            h2d.SetBinContent(ix+1, iy+1, cov_values[ix][iy])
            h2d.SetBinError(ix+1, iy+1, 0)
    return h2d


def same_floats(a, b, rel_tolerance=1E-5):
    diff = abs(a-b)
    return diff/max(1E-100, max(abs(a), abs(b))) < rel_tolerance


def print_tmatrixsparse(matrix, name="matrix"):
    rows = matrix.GetRowIndexArray()
    nrows = matrix.GetNrows()
    cols = matrix.GetColIndexArray()
    ncols = matrix.GetNcols()
    data = matrix.GetMatrixArray()
    for irow in range(nrows):
        sIndex = rows[irow]
        eIndex = rows[irow+1]
        for index in range(sIndex, eIndex):
            icol = cols[index]
            val = data[index]
            print("%s(%d,%d) = %.4e" % (name, irow+matrix.GetRowLwb(), icol+matrix.GetColLwb(), val))


def remove_th1_errors(h):
    for i in range(1, h.GetNbinsX()+1):
        h.SetBinError(i, 0)


def print_th1_bins(h, print_contents=True, print_errors=True, do_oflow=False):
    start_bin = 0 if do_oflow else 1
    end_bin = h.GetNbinsX()+1 if do_oflow else h.GetNbinsX()
    for ix in range(start_bin, end_bin+1):
        print(ix, h.GetBinContent(ix) if print_contents else "", "± %f" % h.GetBinError(ix) if print_errors else "")


def print_th2_bins(h, print_contents=True, print_errors=True, do_oflow=False):
    start_bin_x = 0 if do_oflow else 1
    end_bin_x = h.GetNbinsX()+1 if do_oflow else h.GetNbinsX()
    start_bin_y = 0 if do_oflow else 1
    end_bin_y = h.GetNbinsY()+1 if do_oflow else h.GetNbinsY()
    for ix in range(start_bin_x, end_bin_x+1):
        for iy in range(start_bin_y, end_bin_y+1):
            print(ix, iy, h.GetBinContent(ix. iy) if print_contents else "", "± %f" % h.GetBinError(ix, iy) if print_errors else "")


def get_min_max_bin_contents(hist, include_error=True, ignore_zeros=True, start_bin=None, end_bin=None):
    """Get minimum and maximum bin & contents from hist

    Parameters
    ----------
    hist : ROOT.TH1
        Description
    include_error : bool, optional
        Include error bars in determining minimum & maximum
    ignore_zeros : bool, optional
        If True, then ignore 0s from possible minimum/maximum
    start_bin : None, optional
        Only consider bin numbers >= start_bin.
        If None, start at 1st bin
    end_bin : None, optional
        Only consider bin numbers <= end_bin.
        If None, end at last bin

    Returns
    -------
    float, int, float, int
        Minimum, minimum bin number, maximum, maximum bin number
    """
    min_val = 1E100
    min_bin = 0
    max_val = -1E100
    max_bin = 0
    start_bin = start_bin or 1
    end_bin = end_bin or hist.GetNbinsX()
    for i in range(1, end_bin+1):
        val = hist.GetBinContent(i)
        err = hist.GetBinError(i) if include_error else 0
        val_down = val-err
        val_up = val+err
        if (ignore_zeros and val_down != 0) or not ignore_zeros:
            if (val_down < min_val):
                min_val = val_down
                min_bin = i
        if (ignore_zeros and val_up != 0) or not ignore_zeros:
            if (val_up > max_val):
                max_val = val_up
                max_bin = i
    return min_val, min_bin, max_val, max_bin


def get_min_max_bin_contents_multiple_hists(hists, *args, **kwargs):
    """Do get_min_max_bin_contents across list of hists, getting global min/max & bins

    See get_min_max_bin_contents() for args
    """
    results = [get_min_max_bin_contents(h, *args, **kwargs)
               for h in hists]

    # Get min across all hists, and also the hist bin number corresponding to it
    all_mins = [m[0] for m in results]
    min_all = min(all_mins)
    min_bin = results[all_mins.index(min_all)][1]

    # Get max across all hists, and also the hist bin number corresponding to it
    all_maxs = [m[2] for m in results]
    max_all = max(all_maxs)
    max_bin = results[all_maxs.index(max_all)][3]

    return min_all, min_bin, max_all, max_bin


def get_visible_axis_range(axis):
    # don't use axis.GetXmax(), GetXmin(): they just return the total possible range,
    # not the visible one
    return axis.GetBinLowEdge(axis.GetFirst()), axis.GetFirst(), axis.GetBinLowEdge(axis.GetLast()+1), axis.GetLast()


def get_lumi_str(do_dijet=True, do_zpj=True):
    """Get correct luminosity depending on signal region(s), since dijet prescaled"""
    lumi = "35.9"
    if do_zpj and do_dijet:
        lumi = " #leq35.9"
    elif not do_zpj and do_dijet:
        lumi = " <35.9"
    return lumi

