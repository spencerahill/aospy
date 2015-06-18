import numpy as np

r_e = 6370997. # Earth's mean radius, in meters (m).
c_p = 1003.5 # Specific heat capacity of dry air at constant pressure, in J/K/kg.
c_v = 717. # Specific heat capacity of dry air at constant volume, in J/K/kg.
L_v = 2.5e6 # Latent heat of vaporization of water, in J/kg.
L_f = 3.34e5 # Latent heat of fusion of water, in J/kg.
grav = 9.81 # Acceleration due to gravity, in meters per second squared (m/s^2)
R_d = 287.04 # Dry air gas constant, in J/kg/K
R_v = 461.5 # Water vapor gas constant, in J/kg/K
Omega = 2.*np.pi/(24.*3600.) # Earth's rotation rate.
day2sec = 24.*3600.
sec2day = 1. / day2sec

epsilon = R_d/R_v
kappa = R_d/c_p


