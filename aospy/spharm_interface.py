"""Classes and methods for interfacing with windspharm objects."""
import numpy as np
import spharm
import windspharm

from .constants import r_e

class SpharmInterface(object):
    """Interface between aospy data and spharm and windspharm packages.

    Windspharm has built-in tools for converting to and from common array
    shapes such as (time, p, lat, lon) to the required (lat, lon, records).
    However, the tools do *not* unmask and then re-mask data.  So need to
    write a wrapper that does that part.  Otherwise can just use the
    package's built-in tools.
    """
    def __init__(self, u, v, gridtype='gaussian', rsphere=r_e,
                 legfunc='computed'):
        self.u_raw = u
        self.v_raw = v
        if isinstance(u, np.ma.core.MaskedArray):
            self.mask = u.mask
        else:
            self.mask = False
        self.gridtype = gridtype
        self.rsphere = rsphere
        self.legfunc = legfunc

        self.u, self.grid_info = self.prep_for_spharm(self.u_raw)
        self.v, _ = self.prep_for_spharm(self.v_raw)
        self.n_lat, self.n_lon = self.u.shape[:2]

        self.ws_recover = windspharm.tools.get_recovery(self.grid_info)

        try:
            self.VectorWind = windspharm.standard.VectorWind(
                self.u, self.v, gridtype=self.gridtype
            )
            self.spharmt = self.VectorWind.s
        except (ValueError, AttributeError):
            self.VectorWind = None
            self.spharmt = spharm.Spharmt(
                self.n_lon, self.n_lat, rsphere=self.rsphere,
                gridtype=self.gridtype, legfunc=self.legfunc
            )
        else:
            for func in ('magnitude', 'vrtdiv', 'vorticity', 'divergence',
                         'planetaryvorticity', 'absolutevorticity',
                         'sfvp', 'streamfunction', 'velocitypotential',
                         'helmholtz', 'irrotationalcomponent',
                         'nondivergentcomponent', 'gradient', 'truncate'):
                setattr(self, func, getattr(self.VectorWind, func))

    @staticmethod
    def prep_for_spharm(field):
        filled = np.ma.filled(field, fill_value=0.)[:,:,::-1]
        return windspharm.tools.prep_data(filled, 'tzyx')

    def revert_to_raw(self, field):
        field = self.ws_recover(field)[0]
        field = field[:,:,::-1]
        return np.ma.array(field, mask=self.mask)
