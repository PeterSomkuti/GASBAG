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
maxval = 0.035
#maxval = 5e18
extent = [-maxval, maxval, -maxval, maxval]

# we want to filter out those that are flagged bad in IDP
qual = ((new['physical_retrieval_results/retrieved_chi2_771nm'][:] > 0.001) &
        (new['physical_retrieval_results/retrieved_chi2_771nm'][:] < 2.000))
sounding_ids = new['SoundingGeometry/sounding_id'][:].astype('str')

fig, axarr = plt.subplots(2, 2, figsize=(7, 6))

for i, (idp_key, new_key) in enumerate([#('DOASFluorescence/fluorescence_radiance_757nm_idp',
                                        #f'{retr_type}/retrieved_sif_abs_757nm')]):
                                        #('DOASFluorescence/fluorescence_offset_relative_771nm_idp',
                                        # f'{retr_type}/retrieved_sif_rel_771nm')]):
                                        ('/RetrievalResults/fluorescence_771nm/RetrievedStateVector/state_vector',
                                         f'{retr_type}/retrieved_sif_rel_771nm')]):

    print(i, idp_key, new_key)

    mask = qual & (np.abs(idp[idp_key][:-1,:,5]) < maxval)
    mask = mask & (np.abs(new[new_key][:]) < maxval)
    num_frames = np.sum(mask, axis=1)

    # We want to average over an entire frame here.
    x = np.nanmean(np.ma.masked_array(idp[idp_key][:-1,:,5], mask=~mask), axis=1)[num_frames > 0].compressed()
    y = np.nanmean(np.ma.masked_array(new[new_key][:], mask=~mask), axis=1)[num_frames > 0].compressed()

    #x = np.ma.masked_array(idp[idp_key][:-1,:,5], mask=~mask).compressed()
    #y = np.ma.masked_array(new[new_key][:], mask=~mask).compressed()  #[num_frames > 4].compressed() * 1.8

    print(f"Sum of quality-passed points: {mask.sum()}")

    lreg = sps.linregress(x,y)
    ransac = RANSACRegressor(min_samples=np.floor(0.9*len(x)))
    ransac.fit(x[:, np.newaxis], y[:, np.newaxis])

    print(lreg)

    ax = axarr[0, i]
    ax.set_aspect(1)
    ax.hexbin(x, y, extent=[x.min(), x.max(), y.min(), y.max()],
              linewidths=0, mincnt=1, gridsize=30, cmap='plasma')
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())

    x_ = np.array([-maxval, maxval])

    icept = ransac.estimator_.intercept_[0]
    slope = ransac.estimator_.coef_[0][0]

    labeltext = "f = {:.3g} + ({:.3g} $\pm$ {:.3g}) $\cdot$ x".format(icept, slope, lreg[4])
    ax.plot(x_, x_*slope + icept, 'r--', label=labeltext, lw=1.0)

    axtext = "R = {:.3f}".format(lreg[2])
    ax.text(0.99, 0.05, axtext, ha='right', va='center', transform=ax.transAxes)

    ax.plot(x_, x_, lw=1.0, color='grey', linestyle='dotted', label='1:1')

    if (i == 0):
        ax.set_ylabel("GeoCARB SIF (abs)")
        ax.set_title('757 nm')
    elif (i == 1):
        ax.set_title('771 nm')
    ax.set_xlabel("IDP SIF (abs)")
    ax.legend(fontsize=8)

for i, (idp_key, new_key) in enumerate([('/RetrievalResults/fluorescence_771nm/SpectralParameters/residual_reduced_chi2',
                                         f'{retr_type}/retrieved_chi2_771nm')]):
                                        #('DOASFluorescence/residual_reduced_chi2_fluorescence_771nm_idp',
                                        # f'{retr_type}/retrieved_chi2_771nm')]):

    ax = axarr[1, i]
    ax.hist(idp[idp_key][1:,:][qual].flatten(), range=(0,3), bins=100,
            alpha=0.5, label='IDP')
    ax.hist(new[new_key][:][qual].flatten(), range=(0,3), bins=100,
            alpha=0.5, label='GeoCARB')

    ax.set_xlabel('$\chi^2$')
    ax.legend()

plt.tight_layout(h_pad=2.5, w_pad=2.0)
plt.savefig('idp_comparison_771_2.pdf', bbox_inches='tight')

embed()
