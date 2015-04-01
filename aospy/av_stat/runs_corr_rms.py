import cPickle
import numpy as np
from scipy.stats import pearsonr
from netCDF4 import Dataset
from atmospy import grid_sfc_area

"""Compute rms error and correlation of two runs."""

def weighted_cov(x, y, w):
    """Weighted covariance."""
    return np.sum(w*((x - np.average(x, weights=w)) * 
                     (y - np.average(y, weights=w))))/w.sum()

def weighted_corr(x, y, w):
    """Weighted Pearson correlation coefficient."""
    return (weighted_cov(x, y, w) / 
            np.sqrt(weighted_cov(x, x, w)*weighted_cov(y, y, w)))
            
def weighted_rms(x, y, w):
    """Weighted root-mean-square error."""
    return np.sqrt(np.average((x - y)**2, weights=w))

models = ['am2ml', 'am2']
runs = ['aero-cont','aero-cont']
yrs = ['0061-0080', '1983-1998']
vars = ['precip', 't_surf']
#vars = ('alb_sfc', 'evap', 'high_cld_amt', 'ice_mask', 'IWP', 'low_cld_amt', 'lwdn_sfc_clr', 'lwdn_sfc', 'lwflx', 'LWP', 'lwup_sfc_clr', 'lwup_sfc', 'olr_clr', 'olr', 'prec_conv', 'prec_ls', 'precip', 'ps', 'rh_ref', 'shflx', 'slp', 'snow_conv', 'snow_ls', 'swdn_sfc_clr', 'swdn_sfc', 'swup_sfc_clr', 'swup_sfc', 'swup_toa_clr', 'swup_toa', 't_ref', 't_surf', 'tau_x', 'tau_y', 'tot_cld_amt', 'u_ref', 'v_ref', 'wind', 'WVP')
intvls = ['djf']

pre = '/Users/spencerahill/res/aero_3agcm/data/'
nc = Dataset(pre + 'am2/cont/atmos.1984-1990.ann.nc')
lats = nc.variables['lat'][:]
sfc_area = grid_sfc_area(nc)
sa = sfc_area
#weights = sa
weights = np.abs(np.cos(lats))
#weights = np.ones((90,)) 

rmss = []
corrs = []
p_vals = []

for var in vars:
    print var,
    for intvl in intvls:
        vals = []
        for i, run in enumerate(runs):
            # Load data.
            file = (pre + models[i] + '/' + run + '/' + var + '.av.znl.' + 
                    models[i] + '.' + run + '.' + intvl + '.' + yrs[i] + '.p')
            f = open(file, 'r')
            val = cPickle.load(f)
            # Weight by surface area.
#            val *= sfc_area/sfc_area.sum()
#            vals.append(val.flatten())
            vals.append(val)
        # Compute Pearson correlation coefficient and p value.
        #corr, p_val = pearsonr(vals[0], vals[1])
        corr = weighted_corr(vals[0], vals[1], weights)
        corrs.append(corr)
        #p_vals.append(p_val)
        #print '%0.2f' % corr, p
        # Compute RMS difference.
        rms = weighted_rms(vals[0], vals[1], weights)
        rmss.append(rms)
        print '%0.2f' % corr, '\t', '%0.2f' % rms
corrs, rmss = np.array(corrs), np.array(rmss)
#p_vals = np.array(p_vals)