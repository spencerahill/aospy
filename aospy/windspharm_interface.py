"""Classes and methods for interfacing with windspharm objects."""
import numpy as np
import windspharm


class WindspharmInterface(object):
    """Interface between aospy data and windspharm package.

    Windspharm has built-in tools for converting to and from common array
    shapes such as (time, p, lat, lon) to the required (lat, lon, records).
    However, the tools do *not* unmask and then re-mask data.  So need to
    write a wrapper that does that part.  Otherwise can just use the
    package's built-in tools.
    """
    def __init__(self, u, v, gridtype='regular'):
        self.u_raw = u
        self.v_raw = v
        if isinstance(u, np.ma.core.MaskedArray):
            self.mask = u.mask
        else:
            self.mask = False
        self.gridtype = gridtype

        self.u_prepped, self.grid_info = self.prep_for_windspharm(self.u_raw)
        self.v_prepped, _ = self.prep_for_windspharm(self.v_raw)

        self.ws_recover = windspharm.tools.get_recovery(self.grid_info)

        self.VectorWind = windspharm.standard.VectorWind(
            self.u_prepped, self.v_prepped, gridtype=self.gridtype
        )

        for func in ('magnitude', 'vrtdiv', 'vorticity', 'divergence',
                     'planetaryvorticity', 'absolutevorticity',
                     'sfvp', 'streamfunction', 'velocitypotential',
                     'helmholtz', 'irrotationalcomponent',
                     'nondivergentcomponent', 'gradient', 'truncate'):
            setattr(self, func, getattr(self.VectorWind, func))

    def prep_for_windspharm(self, field):
        filled = np.ma.filled(field, fill_value=0.)[:,:,::-1]
        return windspharm.tools.prep_data(filled, 'tzyx')

    def revert_to_raw(self, field):
        field = self.ws_recover(field)[0]
        field = field[:,:,::-1]
        return np.ma.array(field, mask=self.mask)
