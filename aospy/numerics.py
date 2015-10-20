"""Classes and functions for numerical analysis."""
import numpy as np
import xray

class FiniteDiff(object):
    """For numerical approximations of derivatives using finite differences."""
    def __init__(self, f, geometry='spherical',
                 vector_field=False, wraparound=False):
        """
        Create a `FiniteDiff` object.

        :param f: Data to be finite-differenced.
        :type f: xray.Dataset or xray.DataArray object
        :param str geometry: Geometry of the positions.  Either 'cartesian'
                             or 'spherical'.
        :param vector_field: Whether or not `f` is a component of a vector
                             field.  In some geometries and some directions,
                             operations differ whether the field is scalar or
                             a vector (e.g. north-south on the sphere)
        :param wraparound: Which, if any, axes are wraparound, e.g. longitude
                           on the sphere.
        """
        self.f = f
        self.geometry = geometry
        self.vector_field = vector_field
        self.wraparound = wraparound

    @staticmethod
    def cen_diff(arr, dim, spacing=1):
        """Centered differencing of the DataArray or Dataset.

        :param arr: Data to be center-differenced.
        :type arr: `xray.DataArray` or `xray.Dataset`
        :param str dim: Dimension over which to perform the differencing.
        :param int spacing: How many gridpoints over to use.  Size of resulting
                            array will vary depending on this value.
        """
        if spacing < 1:
            raise ValueError("Centered differencing cannot have spacing < 1")
        left = arr.isel(**{dim: slice(0, -spacing)})
        right = arr.isel(**{dim: slice(spacing, None)})
        return (left.diff(dim, n=1, label='upper') +
                right.diff(dim, n=1, label='lower'))

    @staticmethod
    def fwd_diff1(arr, dim):
        """Forward differencing of the array.  Not its full derivative."""
        return arr.diff(dim, n=1, label='lower')

    @staticmethod
    def bwd_diff1(arr, dim):
        """Backward differencing of the array.  Not its full derivative."""
        return -1*arr.diff(dim, n=1, label='upper')

    @classmethod
    def fwd_diff_deriv(cls, arr, dim):
        """1st order accurate forward differencing approximation of derivative.

        :param arr: Field to take derivative of.
        :param str dim: Values of dimension over which derivative is taken.  If
                        a singleton, assumed to be the uniform spacing.  If an
                        array, assumed to be the values themselves, not the
                        spacing, and must be the same length as `arr`.
        :out: Array containing the df/dx approximation, with length in the 0th
              axis one less than that of the input array.
        """
        return cls.fwd_diff(arr, dim) / cls.fwd_diff(arr[dim], dim)

    @classmethod
    def cen_diff_deriv(cls, arr, dim, order=2, do_edges_one_sided=False):
        """
        Centered differencing approximation of 1st derivative.

        :param arr: Data to be center-differenced.
        :type arr: `xray.DataArray` or `xray.Dataset`
        :param str dim: Dimension over which to compute.
        :param int order: Order of accuracy to use.  Default 2.
        :param do_edges_one_sided: Whether or not to fill in the edge cells
                                   that don't have the needed neighbor cells
                                   for the stencil.  If `True`, use one-sided
                                   differencing with the same order of accuracy
                                   as `order`, and the outputted array is the
                                   same shape as `self.f`.

                                   If `False`, the outputted array has a length
                                   in the computed axis reduced by `order`.
        """
        if order != 2:
            raise NotImplementedError("Centered differencing of df/dx only "
                                      "supported for 2nd order currently")
        deriv = (cls.cen_diff(arr, dim, spacing=1) /
                 cls.cen_diff(arr[dim], dim, spacing=1))
        if do_edges_one_sided:
            left = arr.isel(**{dim: slice(0, 2)})
            deriv_left = (cls.fwd_diff1(left, dim) /
                          cls.fwd_diff1(left[dim], dim))
            right = arr.isel(**{dim: slice(-2, None)})
            deriv_right = (cls.bwd_diff1(right, dim) /
                           cls.bwd_diff1(right[dim], dim))
            deriv = xray.concat([deriv_left, deriv, deriv_right], dim=dim)
        return deriv

    @classmethod
    def upwind_advection(cls, arr, flow, dim, order=1):
        """
        Upwind differencing scheme for advection.

        :param flow: Flow that is advecting the field.
        """
        flow_pos = np.ma.where(flow >= 0., flow, 0)
        flow_neg = np.ma.where(flow < 0., flow, 0)
        return (flow_pos*cls.bwd_diff_deriv(arr, dim, order=order) +
                flow_neg*cls.fwd_diff_deriv(arr, dim, order=order))
