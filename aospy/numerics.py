"""Classes and functions for numerical analysis."""
import numpy as np
import xray


class FiniteDiff(object):
    """For numerical approximations of derivatives using finite differences."""
    def __init__(self, arr, geometry='spherical',
                 vector_field=False, wraparound=False):
        """
        Create a `FiniteDiff` object.

        :param arr: Data to be finite-differenced.
        :type arr: xray.Dataset or xray.DataArray object
        :param str geometry: Geometry of the positions.  Either 'cartesian'
                             or 'spherical'.
        :param vector_field: Whether or not `f` is a component of a vector
                             field.  In some geometries and some directions,
                             operations differ whether the field is scalar or
                             a vector (e.g. north-south on the sphere)
        :param wraparound: Which, if any, axes are wraparound, e.g. longitude
                           on the sphere.
        """
        self.arr = arr
        self.geometry = geometry
        self.vector_field = vector_field
        self.wraparound = wraparound

    @staticmethod
    def fwd_diff1(arr, dim, is_coord=False):
        """Forward differencing of the array.  Not its full derivative.

        A bug in xray version 0.6.1 and prior causes the `DataArray.diff`
        method to not work when applied to a coordinate array.  Therefore,
        a workaround is implemented here and used if the `is_coord` keyword
        argument is True.
        """
        if is_coord:
            arr_diff = arr[dim].diff(dim, n=1, label='lower')
            return xray.DataArray(np.diff(arr[dim]), dims=[dim],
                                  coords=[arr_diff[dim]])
        return arr.diff(dim, n=1, label='lower')

    @staticmethod
    def bwd_diff1(arr, dim, is_coord=False):
        """Backward differencing of the array.  Not its full derivative."""
        if is_coord:
            arr_diff = arr[dim].diff(dim, n=1, label='upper')
            return xray.DataArray(np.diff(arr[dim]), dims=[dim],
                                  coords=[arr_diff[dim]])
        return arr.diff(dim, n=1, label='upper')

    @classmethod
    def cen_diff(cls, arr, dim, spacing=1, is_coord=False,
                 do_edges_one_sided=False):
        """Centered differencing of the DataArray or Dataset.

        :param arr: Data to be center-differenced.
        :type arr: `xray.DataArray` or `xray.Dataset`
        :param str dim: Dimension over which to perform the differencing.
        :param int spacing: How many gridpoints over to use.  Size of resulting
                            array depends on this value.
        :param do_edges_one_sided: Whether or not to fill in the edge cells
                                   that don't have the needed neighbor cells
                                   for the stencil.  If `True`, use one-sided
                                   differencing with the same order of accuracy
                                   as `order`, and the outputted array is the
                                   same shape as `arr`.

                                   If `False`, the outputted array has a length
                                   in the computed axis reduced by `order`.
        """
        if spacing < 1:
            raise ValueError("Centered differencing cannot have spacing < 1")
        left = arr.isel(**{dim: slice(0, -spacing)})
        right = arr.isel(**{dim: slice(spacing, None)})
        # Centered differencing = sum of intermediate forward differences
        diff = (cls.fwd_diff1(right, dim, is_coord=is_coord) +
                cls.bwd_diff1(left, dim, is_coord=is_coord))
        if do_edges_one_sided:
            left = arr.isel(**{dim: slice(0, 2)})
            right = arr.isel(**{dim: slice(-2, None)})
            diff_left = cls.fwd_diff1(left, dim, is_coord=is_coord)
            diff_right = cls.bwd_diff1(right, dim, is_coord=is_coord)
            diff = xray.concat([diff_left, diff, diff_right], dim=dim)
        return diff

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
        return (cls.fwd_diff1(arr, dim, is_coord=False) /
                cls.fwd_diff1(arr[dim], dim, is_coord=True))

    @classmethod
    def bwd_diff_deriv(cls, arr, dim):
        """1st order accurate backward differencing approx of derivative.

        :param arr: Field to take derivative of.
        :param str dim: Values of dimension over which derivative is taken.  If
                        a singleton, assumed to be the uniform spacing.  If an
                        array, assumed to be the values themselves, not the
                        spacing, and must be the same length as `arr`.
        :out: Array containing the df/dx approximation, with length in the 0th
              axis one less than that of the input array.
        """
        return (cls.bwd_diff1(arr, dim, is_coord=False) /
                cls.bwd_diff1(arr[dim], dim, is_coord=True))

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
                                   same shape as `arr`.

                                   If `False`, the outputted array has a length
                                   in the computed axis reduced by `order`.
        """
        if order != 2:
            raise NotImplementedError("Centered differencing of df/dx only "
                                      "supported for 2nd order currently")
        return (cls.cen_diff(arr, dim, spacing=1, is_coord=False,
                             do_edges_one_sided=do_edges_one_sided) /
                cls.cen_diff(arr[dim], dim, spacing=1, is_coord=True,
                             do_edges_one_sided=do_edges_one_sided))

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
