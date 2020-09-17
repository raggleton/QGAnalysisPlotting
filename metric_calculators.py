"""
Various functions to calculate metrics

Different methods depending on package, input type, etc

Yes it's a bit of a mess
"""


from __future__ import print_function, division

import numpy as np
import warnings


def get_hist_bin_widths(hist):
    """Get bin widths as numpy array

    Parameters
    ----------
    hist : uproot.TH1

    Returns
    -------
    numpy.array
    """
    return hist.edges[1:] - hist.edges[:-1]


def uproot_th1_to_arrays(hist):
    """Convert histogram bins to arrays

    Assumes `hist` is an uproot.TH1 object, not a ROOT.TH1

    Note that this returns bin *areas* and not heights: assume bin contents
    already divided by bin width.

    Note that errors here are multiplied by the bin width (since they are
    assumed to have been created by originally dividing by the bin width)
    """
    bin_widths = get_hist_bin_widths(hist)
    bin_areas = hist.values * bin_widths
    bin_centers = hist.edges[:-1] + (0.5*bin_widths)
    bin_errors = np.sqrt(hist.variances) * bin_widths
    return bin_areas, bin_widths, bin_centers, bin_errors

try:
    # jax for differentiation to do errors on uproot.TH1 NOT ROOT.TH1
    # NB would like to replace with uncertainties package,
    # but covariance matrix handling tricky
    import jax.numpy as jnp
    from jax import grad, jit

    # from jax.config import config
    # to debug NaNs - turn off if not debugging as slow
    # config.update("jax_debug_nans", True)

    def calc_mean_jax(bin_areas, bin_centers):
        """Calculate mean of hist from value arrays.

        Assumes inputs are numpy

        Must use jnp.X functions for e.g. sum(), square(), to ensure jax can differentiate it
        """
        return jnp.sum(bin_areas * bin_centers) / jnp.sum(bin_areas)


    mean_differential_jax = jit(grad(calc_mean_jax, argnums=0))


    def calc_mean_uncorrelated_error_jax(bin_areas, bin_centers, bin_errors):
        """Calculate error on mean, assuming uncorrelated errors.

        Uses propagation of uncertainty bia partial differentials,
        calculated automatically using jax.
        """
        # differential wrt bin_areas
        diffs = mean_differential_jax(bin_areas, bin_centers)
        err_sq = jnp.sum(jnp.square((diffs * bin_errors)))
        return jnp.sqrt(err_sq)


    def calc_mean_cov_matrix_jax(bin_areas, bin_centers, error_matrix):
        """Get error on mean, assuming covariance matrix error_matrix"""
        diffs = mean_differential_jax(bin_areas, bin_centers)
        # sum_sq = diffs.T @ error_matrix @ diffs
        sum_sq = jnp.matmul(diffs.T, jnp.matmul(error_matrix, diffs))
        return jnp.sqrt(sum_sq)


    def calc_hist_mean_and_uncorrelated_error_jax(hist):
        """Calculate from hist both mean and its error,
        assuming uncorrelated uncertainties

        Parameters
        ----------
        hist : uproot.TH1

        Returns
        -------
        float, float
        """
        areas, widths, centers, errors = uproot_th1_to_arrays(hist)
        mean = calc_mean_jax(areas, centers)
        err = calc_mean_uncorrelated_error_jax(areas, centers, errors)
        return float(mean), float(err)


    def calc_hist_mean_and_correlated_error_jax(hist, ematrix):
        """Calculate from hist both mean and its error,
        uncertainties in the form of a covariance matrix

        Parameters
        ----------
        hist : uproot.TH1
            Description
        ematrix : TYPE
            Description
        hist : uproot.TH1

        Returns
        -------
        float, float

        """
        areas, widths, centers, errors = uproot_th1_to_arrays(hist)
        mean = calc_mean_jax(areas, centers)
        err = calc_mean_cov_matrix_jax(areas, centers, ematrix)
        return float(mean), float(err)

    # --------------------------------
    # Functions to calculate RMS
    # --------------------------------

    def calc_rms_jax(bin_areas, bin_centers):
        """Calculate RMS of hist from value arrays.

        Must use jnp.X functions for e.g. sum(), square(), to ensure jax can differentiate it
        """
        mean = calc_mean_jax(bin_areas, bin_centers)
        sum_sq = jnp.sum(jnp.square((bin_areas * bin_centers) - mean))
        return jnp.sqrt(sum_sq / jnp.sum(bin_areas))


    rms_differential_jax = jit(grad(calc_rms_jax, argnums=0))


    def calc_rms_uncorrelated_error_jax(bin_areas, bin_centers, bin_errors):
        """Calculate error on RMS, assuming uncorrelated errors.

        Uses propagation of uncertainty bia partial differentials,
        calculated automatically using jax.
        """
        # differential wrt bin_areas
        diffs = rms_differential_jax(bin_areas, bin_centers)
        err_sq = jnp.sum(jnp.square((diffs * bin_errors)))
        return jnp.sqrt(err_sq)


    def calc_rms_correlated_error_jax(bin_areas, bin_centers, error_matrix):
        """Get error on rms, assuming covariance matrix error_matrix"""
        diffs = rms_differential_jax(bin_areas, bin_centers)
        # sum_sq = diffs @ error_matrix @ diffs
        sum_sq = jnp.matmul(diffs, jnp.matmul(error_matrix, diffs))
        return jnp.sqrt(sum_sq)


    def calc_hist_rms_and_uncorrelated_error_jax(hist):
        """Calculate from hist both RMS and its error,
        assuming uncorrelated uncertainties

        Parameters
        ----------
        hist : uproot.TH1

        Returns
        -------
        float, float
        """
        areas, widths, centers, errors = uproot_th1_to_arrays(hist)
        rms = calc_rms_jax(areas, centers)
        err = calc_rms_uncorrelated_error_jax(areas, centers, errors)
        return float(rms), float(err)


    def calc_hist_rms_and_correlated_error_jax(hist, ematrix):
        """Summary

        Parameters
        ----------
        hist : TYPE
            Description
        ematrix : TYPE
            Description

        Returns
        -------
        TYPE
            Description
        """
        areas, widths, centers, errors = uproot_th1_to_arrays(hist)
        rms = calc_rms_jax(areas, centers)
        err = calc_rms_correlated_error_jax(areas, centers, ematrix)
        return float(rms), float(err)

    # --------------------------------
    # Functions to calculate delta
    # --------------------------------
    def calc_delta_jax(areas_a, areas_b):
        # do I need bin areas or densities?
        # I guess since by definition sum(area_a) = 1, areas are needed?!
        integrand = jnp.true_divide(jnp.square(areas_a - areas_b), areas_a + areas_b)
        # nan_to_num important as divide gives nans if both 0
        delta = 0.5 * jnp.sum(jnp.nan_to_num(integrand))
        return delta


    delta_differential_jax = jit(grad(calc_delta_jax, argnums=[0, 1]))


    def calc_delta_uncorrelated_error_jax(areas_a, errors_a, areas_b, errors_b):
        pass


    def calc_delta_correlated_error_jax(areas_a, ematrix_a, areas_b, errors_b):
        diffs_a, diffs_b = delta_differential_jax(areas_a, areas_b)
        # need to do nan_to_num since the differential can return nan...
        # not sure how to fix "properly" though
        diffs_a = jnp.nan_to_num(diffs_a)
        diffs_b = jnp.nan_to_num(diffs_b)
        # for the total, we need to do
        # diffs_a * ematrix_a * diffs_a + diffs_b*errors_b*diffs_b,
        # since the errors on a and b have no connections, we can get away with this.
        # err_a_sq = diffs_a.T @ ematrix_a @ diffs_a
        err_a_sq = jnp.matmul(diffs_a.T, jnp.matmul(ematrix_a, diffs_a))
        err_b_sq = jnp.sum(jnp.square((diffs_b * errors_b)))
        return jnp.sqrt(err_a_sq + err_b_sq)


    def calc_hist_delta_and_error_jax(hist_a, ematrix_a, hist_b):
        """Calculate delta between hists, along with its error

        Defined as 0.5 * integral[ (a - b)^2 / (a+b) ]
        """
        areas_a, widths_a, centers_a, errors_a = uproot_th1_to_arrays(hist_a)
        areas_b, widths_b, centers_b, errors_b = uproot_th1_to_arrays(hist_b)
        delta = calc_delta_jax(areas_a, areas_b)
        err = calc_delta_correlated_error_jax(areas_a, ematrix_a, areas_b, errors_b)
        return float(delta), float(err)

except ImportError:
    warnings.warn("jax package is not available: cannot use *_jax() functions")


try:
    import yoda
    # uncertainties package for RIVET/YODA plots
    # since NAF doesn't have any other auto-diff packages,
    # and I cba to write out all the differentials manually
    import uncertainties
    from uncertainties import unumpy as unp


    def yoda_hist_to_values(hist):
        """Get bin centers, widths, heights & errors as individual numpy arrays

        FIXME: make consistent wih the other method

        Parameters
        ----------
        hist : yoda.Histo1D or yoda.Scatter2D

        Returns
        -------
        np.array, np.array, np.array, np.array
        """
        # return bin_areas, bin_widths, bin_centers, bin_errors

        if isinstance(hist, yoda.core.Histo1D):
            centers = np.array([b.xMid for b in hist.bins])
            widths = np.array([b.xWidth for b in hist.bins])
            areas = np.array([b.area for b in hist.bins])
            errors = np.array([b.areaErr for b in hist.bins])
            # return centers, widths, heights, errors
            return areas, widths, centers, errors
        elif isinstance(hist, yoda.core.Scatter2D):
            centers = np.array([p.x for p in hist.points])
            widths = np.array([p.xErrs.plus + p.xErrs.minus for p in hist.points])
            heights = np.array([p.y for p in hist.points])
            areas = heights * widths
            errors = np.array([p.yErrs.plus for p in hist.points])  # assumes symmetric error bars!
            errors = errors * widths
            # return centers, widths, heights, errors
            return areas, widths, centers, errors
        else:
            raise TypeError('hist_to_values only accepts Histo1D or Scatter2D objects, not %s' % (type(hist)))


    def hist_values_to_uarray(bin_areas, bin_centers, bin_errors):
        """Convert numpy arrays of areas and centers to uncertainties.unumpy.uarray objects,
        which store value +- error and can be used in *_ucert() methods

        Inputs comes from yoda_hist_to_values()

        Parameters
        ----------
        bin_areas : numpy.array
        bin_centers : numpy.array
        bin_errors : numpy.array
        """
        areas = unp.uarray(bin_areas, bin_errors)
        centers = unp.uarray(bin_centers, np.zeros_like(bin_centers))
        return areas, centers


    def calc_mean_ucert(bin_areas, bin_centers):
        """Calculate mean of hist from value arrays

        Input arrays are uncertainties.unumpy.uarray
        """
        return (bin_areas * bin_centers).sum() / bin_areas.sum()


    def calc_rms_ucert(bin_areas, bin_centers):
        """Calculate RMS of hist from value arrays

        Input arrays are uncertainties.unumpy.uarray


        Must use uncertainties functions for e.g. square()
        """
        mean = calc_mean_ucert(bin_areas, bin_centers)
        sum_sq = (((bin_areas * bin_centers) - mean)**2).sum()
        return uncertainties.umath_core.sqrt(sum_sq / bin_areas.sum())


    def calc_delta_ucert(bin_areas_a, bin_areas_b):
        """Calculate delta of hist from value arrays

        Input arrays are uncertainties.unumpy.uarray
        """
        this_bin_areas_a = bin_areas_a
        if bin_areas_a.sum().nominal_value != 1:
            print("bin_areas_a not normalised to 1, will normalise sum to 1")
            this_bin_areas_a = bin_areas_a / bin_areas_a.sum().nominal_value
        this_bin_areas_b = bin_areas_b
        if bin_areas_b.sum().nominal_value != 1:
            print("bin_areas_b not normalised to 1, will normalise sum to 1")
            this_bin_areas_b = bin_areas_b / bin_areas_b.sum().nominal_value
        sum_bins = this_bin_areas_a + this_bin_areas_b
        mask = unp.nominal_values(sum_bins) != 0
        diff_sq = (this_bin_areas_a - this_bin_areas_b)**2
        integrand = diff_sq[mask] / sum_bins[mask]
        delta = 0.5 * integrand.sum()
        return delta

except ImportError:
    warnings.warn("uncertainties and/or yoda package is not available: cannot use *_ucert() functions")
