import h5py
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats as sps

from sklearn.linear_model import RANSACRegressor


from IPython import embed
import sys

def my_mode(x):
    return sps.mode(x)[0][0]

new_fname = sys.argv[1]
idp_fname = sys.argv[2]

idp = h5py.File(idp_fname, 'r')
new = h5py.File(new_fname, 'r')

retr_type = "physical_retrieval_results"

# extent
maxval = 0.03
#maxval = 3.5e18
extent = [-maxval, maxval, -maxval, maxval]

# we want to filter out those that are flagged bad in IDP
qual = ((new['physical_retrieval_results/retrieved_chi2_771nm'][:] > 0.001) &
        (new['physical_retrieval_results/retrieved_chi2_771nm'][:] < 2.000))
sounding_ids = new['SoundingGeometry/sounding_id'][:].astype('str')

fig, axarr = plt.subplots(2, 2, figsize=(7, 6))

for i, (idp_key, new_key) in enumerate([#('DOASFluorescence/fluorescence_radiance_771nm_idp',
                                        # f'{retr_type}/retrieved_sif_abs_771nm')]):
                                        ('DOASFluorescence/fluorescence_offset_relative_771nm_idp',
                                         f'{retr_type}/retrieved_sif_rel_771nm')]):

    print(i, idp_key, new_key)

    mask = qual & (np.abs(idp[idp_key][:]) < maxval)
    mask = mask & (np.abs(new[new_key][:]) < maxval)

    print(f"Sum of quality-passed points: {mask.sum()}")

    fp = np.array([int(x[-1]) for x in sounding_ids[mask]])

    lreg = sps.linregress(idp[idp_key][:][mask],
                          new[new_key][:][mask])

    ransac = RANSACRegressor(min_samples=np.floor(0.9*mask.sum()))
    ransac.fit(idp[idp_key][:][mask][:, np.newaxis], new[new_key][:][mask][:, np.newaxis])

    for this_fp in range(1,9):
        lreg2 = sps.linregress(idp[idp_key][:][mask][this_fp == fp],
                               new[new_key][:][mask][this_fp == fp])
        print(this_fp, lreg2)

    print(lreg)

    ax = axarr[0, i]
    ax.set_aspect(1)
    ax.hexbin(idp[idp_key][:][mask],
              new[new_key][:][mask],
              C=idp['DOASFluorescence/continuum_level_radiance_771nm_idp'][:][mask],
              extent=extent, linewidths=0, gridsize=100, mincnt=1, cmap='plasma')

    x = np.array([-maxval, maxval])

    icept = ransac.estimator_.intercept_[0]
    slope = ransac.estimator_.coef_[0][0]

    labeltext = "f = {:.3g} + ({:.3g} $\pm$ {:.3g}) $\cdot$ x".format(icept, slope, lreg[4])
    ax.plot(x, x*slope + icept, 'r--', label=labeltext, lw=1.0)

    axtext = "R = {:.3f}".format(lreg[2])
    ax.text(0.99, 0.05, axtext, ha='right', va='center', transform=ax.transAxes)

    ax.plot(x, x, lw=1.0, color='grey', linestyle='dotted', label='1:1')

    if (i == 0):
        ax.set_ylabel("GeoCARB SIF (abs)")
        ax.set_title('757 nm')
    elif (i == 1):
        ax.set_title('771 nm')
    ax.set_xlabel("IDP SIF (abs)")
    ax.legend(fontsize=8)

for i, (idp_key, new_key) in enumerate([#('DOASFluorescence/residual_reduced_chi2_fluorescence_771nm_idp',
                                        # f'{retr_type}/retrieved_chi2_771nm')]):
                                        ('DOASFluorescence/residual_reduced_chi2_fluorescence_771nm_idp',
                                         f'{retr_type}/retrieved_chi2_771nm')]):

    ax = axarr[1, i]
    ax.hist(idp[idp_key][:][qual].flatten(), range=(0,3), bins=100,
            alpha=0.5, label='IDP')
    ax.hist(new[new_key][:][qual].flatten(), range=(0,3), bins=100,
            alpha=0.5, label='GeoCARB')

    ax.set_xlabel('$\chi^2$')
    ax.legend()

plt.tight_layout(h_pad=2.5, w_pad=2.0)
plt.savefig('idp_comparison.pdf', bbox_inches='tight')
