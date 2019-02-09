import h5py
import numpy as np
from matplotlib import pyplot as plt
import sys

SH_H2O_CONV = 1.607524279389496


f = h5py.File(sys.argv[1], 'r')
atmos = np.genfromtxt(sys.argv[2], skip_header=1)


p_levels_met = f['ECMWF/vector_pressure_levels_ecmwf'][:,0,:]
sh_met = f['ECMWF/specific_humidity_profile_ecmwf'][:,0,:]

sh_met = sh_met / (1.0 - sh_met) * SH_H2O_CONV

sh_on_atmos = np.zeros((sh_met.shape[0], atmos.shape[0]))

for i in range(sh_on_atmos.shape[0]):
    sh_on_atmos[i] = np.interp(atmos[:,0], p_levels_met[i], sh_met[i])

cov = np.clip(np.cov(sh_on_atmos.T), 0, 1e10)

#cov[cov < 1e-10] = 0.0

# Set diagonals to at least some small value
for i in range(cov.shape[0]):
    if (cov[i,i] < 1e-10):
        cov[i,i] = 1e-6

cov *= 1.0

np.savetxt('h2o_covariance.dat', cov, fmt='%.5e')

