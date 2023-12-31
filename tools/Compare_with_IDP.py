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


sounding_ids = new['SoundingGeometry/sounding_id'][:].astype('str')

fig, axarr = plt.subplots(2, 2, figsize=(7, 6))

for i, (idp_key, new_key) in enumerate([#('DOASFluorescence/fluorescence_radiance_757nm_idp',
                                        #f'{retr_type}/retrieved_sif_abs_757nm')]):
                                        #('DOASFluorescence/fluorescence_offset_relative_771nm_idp',
                                        # f'{retr_type}/retrieved_sif_rel_771nm')]):
        ('/RetrievalResults/fluorescence_757nm/RetrievedStateVector/state_vector',
         f'{retr_type}/retrieved_sif_rel_757nm'),
        ('/RetrievalResults/fluorescence_771nm/RetrievedStateVector/state_vector',
         f'{retr_type}/retrieved_sif_rel_771nm')]):


    if (i==0):
        win = '757nm'
    elif (i==1):
        win = '771nm'

    print(i, idp_key, new_key)

    # we want to filter out those that are flagged bad in IDP
    qual = (np.isfinite(new[f'{retr_type}/retrieved_chi2_{win}'][:]))
    qual = (qual &
            (new[f'physical_retrieval_results/retrieved_chi2_{win}'][:] < 3.00) &
            (new[f'physical_retrieval_results/retrieved_num_iterations_{win}'][:] < 4))

    if (i == 0):
        mask = qual & (np.abs(idp[idp_key][:-1,:,4]) < maxval)
        x = idp[idp_key][:-1,:,4][mask]
    elif (i==1):
        mask = qual & (np.abs(idp[idp_key][:-1,:,5]) < maxval)
        x = idp[idp_key][:-1,:,5][mask]

    #mask = mask & (np.abs(new[new_key][:]) < maxval)
    #num_frames = np.sum(mask, axis=1)

    # We want to average over an entire frame here.
    y = new[f'{retr_type}/retrieved_sif_abs_{win}'][:] / (
        idp[f'RetrievalResults/fluorescence_{win}/Ancillary/continuumLevelRadiance'][:-1,:] - new[f'{retr_type}/retrieved_sif_abs_{win}'][:])
    #y = np.nanmean(np.ma.masked_array(y, mask=~mask), axis=1).compressed()
    y = y[mask]
    #y = np.nanmean(np.ma.masked_array(new[f'{retr_type}/retrieved_sif_rel_771nm'][:], mask=~mask), axis=1)[num_frames > 3].compressed()

    #x = np.ma.masked_array(idp[idp_key][:-1,:,5], mask=~mask).compressed()
    #y = np.ma.masked_array(new[new_key][:], mask=~mask).compressed()  #[num_frames > 4].compressed() * 1.8

    print(f"Sum of quality-passed points: {mask.sum()}")

    lreg = sps.linregress(x,y)
    ransac = RANSACRegressor(min_samples=np.floor(0.9*len(x)))
    ransac.fit(x[:, np.newaxis], y[:, np.newaxis])

    print(lreg)

    ax = axarr[0, i]
    ax.set_aspect(1)
    ex_min = min(x.min(), y.min())
    ex_max = max(x.max(), y.max())
    ax.hexbin(x, y, extent=[ex_min, ex_max, ex_min, ex_max], #extent=[x.min(), x.max(), y.min(), y.max()],
              linewidths=0, mincnt=2, gridsize=50, cmap='plasma')
    ax.set_xlim(ex_min, ex_max)
    ax.set_ylim(ex_min, ex_max)

    x_ = np.array([ex_min, ex_max])

    icept = ransac.estimator_.intercept_[0]
    slope = ransac.estimator_.coef_[0][0]

    labeltext = "f = {:.3g} + ({:.3g} $\pm$ {:.3g}) $\cdot$ x".format(icept, slope, lreg[4])
    ax.plot(x_, x_*slope + icept, 'r--', label=labeltext, lw=1.0)

    axtext = "R = {:.3f}".format(lreg[2])
    ax.text(0.99, 0.05, axtext, ha='right', va='center', transform=ax.transAxes)

    ax.plot(x_, x_, lw=1.0, color='grey', linestyle='dotted', label='1:1')

    ax.grid()

    if (i == 0):
        ax.set_ylabel("GeoCARB SIF (rel)")
        ax.set_title('757 nm')
    elif (i == 1):
        ax.set_title('771 nm')
    ax.set_xlabel("IDP SIF (rel)")
    ax.legend(fontsize=6, loc='upper left')

for i, (idp_key, new_key) in enumerate([('/RetrievalResults/fluorescence_757nm/SpectralParameters/residual_reduced_chi2',
                                         f'{retr_type}/retrieved_chi2_757nm'),
                                        ('/RetrievalResults/fluorescence_771nm/SpectralParameters/residual_reduced_chi2',
                                         f'{retr_type}/retrieved_chi2_771nm')]):
                                        #('DOASFluorescence/residual_reduced_chi2_fluorescence_771nm_idp',
                                        # f'{retr_type}/retrieved_chi2_771nm')]):

    if (i==0):
        win = '757nm'
    elif (i==1):
        win = '771nm'

    qual = ((new[f'{retr_type}/retrieved_chi2_{win}'][:] > 0.001) &
            (new[f'physical_retrieval_results/retrieved_chi2_{win}'][:] < 3.500))

    ax = axarr[1, i]
    ax.hist(idp[idp_key][1:,:][qual].flatten(), range=(0,3), bins=100,
            alpha=0.5, label='IDP')
    ax.hist(new[new_key][:][qual].flatten(), range=(0,3), bins=100,
            alpha=0.5, label='GeoCARB')

    if (i==0):
        ax.set_ylabel('Fit quality')

    ax.set_xlabel('$\chi^2$')
    ax.legend()

plt.tight_layout(h_pad=2.5, w_pad=2.0)
plt.savefig('idp_comparison_both.pdf', bbox_inches='tight')

embed()
