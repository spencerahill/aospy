import numpy as np
import numpy.ma as ma
from scipy.stats import scoreatpercentile
import matplotlib.cm as cm
import matplotlib.colors as colors

def subcolorbar(xmin, xmax, cbar):
    """
    Adapted by Spencer Hill from script created by Chloe Lewis. Returns the part
    of cbar between xmin, xmax, scaled to 0,1.
    """
    assert xmin < xmax
    assert xmax <= 1
    cd = cbar._segmentdata.copy()
    rgbmin = dict(zip(('red', 'green', 'blue', 'alpha'), cbar(xmin)))
    rgbmax = dict(zip(('red', 'green', 'blue', 'alpha'), cbar(xmax)))
    for k in cd:
        tmp = [x for x in cd[k] if x[0] >= xmin and x[0] <= xmax]
        if tmp == [] or tmp[0][0] > xmin:
            tmp = [(xmin, rgbmin[k], rgbmin[k])] + tmp
        if tmp == [] or tmp[-1][0] < xmax:
            tmp = tmp + [ (xmax,rgbmax[k], rgbmax[k])]
        # Now scale all this to (0,1).
        square = zip(*tmp)
        xbreaks = [(x - xmin)/(xmax-xmin) for x in square[0]]
        square[0] = xbreaks
        tmp = zip(*square)
        cd[k] = tmp
    return colors.LinearSegmentedColormap('local', cd, N=256)

def trunc_col_map(data, thresh, n_cntr, col_map):
    """
    Create color map centered on zero and extending to given percentiles.
    """
    data = ma.compressed(data)
    upper = np.ceil(scoreatpercentile(data.flatten(), 100.-thresh))
    lower = np.floor(scoreatpercentile(data.flatten(), thresh))
    low_norm = 0.5 - np.abs(lower/upper)*0.5
    up_norm = 0.5 + np.abs(upper/lower)*0.5
    cntr_lvls = np.linspace(lower, upper, n_cntr)
    # Color map to use.
    cbar_cmd = 'cbar = cm.' + col_map
    exec(cbar_cmd)
    # Make color map, based on whether more data is positive or negative.
    if np.abs(upper/lower) > 1:
        cmap = subcolorbar(low_norm, 1, cbar)
    else:
        cmap = subcolorbar(0, up_norm, cbar)
    return cmap, cntr_lvls