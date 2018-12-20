import h5py
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats as sps

from IPython import embed
import sys

def my_mode(x):
    return sps.mode(x)[0][0]

new_fname = sys.argv[1]
idp_fname = sys.argv[2]

idp = h5py.File(idp_fname, 'r')
new = h5py.File(new_fname, 'r')

# extent
maxval = 0.085
extent = [-maxval, maxval, -maxval, maxval]

# we want to filter out those that are flagged bad in IDP
qual = idp['DOASFluorescence/fluorescence_qual_flag_idp'][:] == 0

sounding_ids = new['SoundingGeometry/sounding_id'][:].astype('str')

fig, axarr = plt.subplots(2, 2, figsize=(7, 6))

for i, (idp_key, new_key) in enumerate([('DOASFluorescence/fluorescence_offset_relative_757nm_idp',
                                         'linear_fluorescence_results/retrieved_sif_rel_757nm'),
                                        ('DOASFluorescence/fluorescence_offset_relative_771nm_idp',
                                         'linear_fluorescence_results/retrieved_sif_rel_771nm')]):

    print(i, idp_key, new_key)

    mask = qual & (np.abs(idp[idp_key][:] < maxval))
    mask = mask & (np.abs(new[new_key][:] < maxval))

    fp = np.array([int(x[-1]) for x in sounding_ids[mask]])

    lreg = sps.linregress(idp[idp_key][:][mask],
                          new[new_key][:][mask])

    ax = axarr[0, i]
    ax.set_aspect(1)
    ax.hexbin(idp[idp_key][:][mask],
              new[new_key][:][mask],
              #C=fp, reduce_C_function=my_mode,
              extent=extent, linewidths=0, gridsize=100, mincnt=1, cmap='plasma')

    x = np.array([-maxval, maxval])
    labeltext = "f = {:.3f} + ({:.3f} $\pm$ {:.3f}) $\cdot$ x".format(lreg[1], lreg[0], lreg[4])
    ax.plot(x, x*lreg[0] + lreg[1], 'r--', label=labeltext, lw=1.0)

    axtext = "R = {:.3f}".format(lreg[2])
    ax.text(0.99, 0.05, axtext, ha='right', va='center', transform=ax.transAxes)

    ax.plot(x, x, lw=1.0, color='grey', linestyle='dotted', label='1:1')

    if (i == 0):
        ax.set_ylabel("GeoCARB SIF (relative)")
        ax.set_title('757 nm')
    elif (i == 1):
        ax.set_title('771 nm')
    ax.set_xlabel("IDP SIF (relative)")
    ax.legend(fontsize=8)

for i, (idp_key, new_key) in enumerate([('DOASFluorescence/residual_reduced_chi2_fluorescence_757nm_idp',
                                         'linear_fluorescence_results/retrieved_chi2_757nm'),
                                        ('DOASFluorescence/residual_reduced_chi2_fluorescence_771nm_idp',
                                         'linear_fluorescence_results/retrieved_chi2_771nm')]):

    ax = axarr[1, i]
    ax.hist(idp[idp_key][:].flatten(), range=(0,3), bins=100,
            alpha=0.5, label='IDP')
    ax.hist(new[new_key][:].flatten(), range=(0,3), bins=100,
            alpha=0.5, label='GeoCARB')

    ax.set_xlabel('$\chi^2$')
    ax.legend()

plt.tight_layout(h_pad=2.5, w_pad=2.0)
plt.savefig('idp_comparison.pdf', bbox_inches='tight')
