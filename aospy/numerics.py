"""Classes and functions for numerical analysis."""
import numpy as np


class FiniteDiff(object):
    """For numerical approximations of derivatives using finite differences."""
    def __init__(self, field, positions=False, geometry='spherical',
                 vector_field=False, wraparound=False):
        """
        Create a `FiniteDiff` object.

        :param field: Field of data to be finite-differenced.
        :param positions: Array of positions (e.g. in space or time) used
                          to determine grid spacing.
        :param geometry: Geometry of the positions.  Options are 'cartesian'
                         or 'spherical'.
        :param vector_field: Whether or not `field` is a component of a vector
                             field.  In some geometries and some directions,
                             operations differ whether the field is scalar or
                             a vector (e.g. north-south on the sphere)
        :param wraparound: Which, if any, axes are wraparound, e.g. longitude
                           on the sphere.
        """
        self.field = field
        if positions is False:
            self.positions = np.arange(len(field)).reshape(field.shape)
        elif isinstance(positions, (float, int)):
            self.positions = np.arange(len(field)).reshape(field.shape)*positions
        else:
            self.positions = positions
        self.geometry = geometry
        self.vector_field = vector_field
        self.wraparound = wraparound

    def cen_diff(self, array, axis=0, spacing=2):
        """Centered differencing of the array.  NOT its full derivative."""
        if spacing < 2:
            raise ValueError("Centered differencing cannot have spacing < 2")
        if axis:
            raise NotImplementedError("Only 0th axis supported for now.")
        return array[spacing:] - array[:-spacing]

    def fwd_diff(self, array, axis=0, spacing=1):
        """Forward differencing of the array.  NOT its full derivative."""
        if spacing < 1:
            raise ValueError("Forward and backward differencing cannot have "
                             "spacing < 1")
        if axis:
            raise NotImplementedError("Only 0th axis supported for now.")
        return array[spacing:] - array[:-spacing]

    def bwd_diff(self, array, axis=0, spacing=1):
        """Backward differencing of the array.  Not its full derivative."""
        return self.fwd_diff(array[::-1])[::-1]

    def cen_diff_deriv(self, axis=0, order=2, fill_edges=False):
        """
        Centered differencing approximation of 1st derivative.

        :param axis: Axis number over which to compute.  (Currently only
                     axis=0 is supported.)
        :param order: Order of accuracy to use.  Default 2.
        :param fill_edges: Whether or not to fill in the edge cells that don't
                           have the needed neighbor cells for the stencil.  If
                           `True`, use one-sided differencing with the same
                           order of accuracy as `order`, and the outputted
                           array is the same shape as `self.field`.

                           If `False`, the outputted array has a length in
                           the computed axis reduced by `order`.
        """
        if order == 2:
            return (self.cen_diff(self.field, axis=axis, spacing=2) /
                    self.cen_diff(self.positions, axis=axis, spacing=2))
        else:
            raise NotImplementedError("Centered differencing of df/dx only "
                                      "supported for 2nd order currently")

    def upwind_advection(self, flow, axis=0, order=1):
        """
        Upwind differencing scheme for advection.

        :param flow: Flow that is advecting the field.
        """
        flow_pos = np.ma.where(flow >= 0., flow, 0)
        flow_neg = np.ma.where(flow < 0., flow, 0)
        return (flow_pos*self.bwd_diff_deriv(axis=axis, order=order) +
                flow_neg*self.fwd_diff_deriv(axis=axis, order=order))
