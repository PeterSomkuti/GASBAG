from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times', 'Times New Roman']
rcParams['mathtext.fontset'] = 'stix'


import numpy as np
from matplotlib import pyplot as plt
import scipy as sp
from scipy import stats as sps
import h5py
import sys

# Grab the file names from command line arguments
if len(sys.argv) != 3:
    print("Usage: python Compare_with_IDP_v2.py new.h5 old.h5")
    sys.exit()

new_fname = sys.argv[1]
idp_fname = sys.argv[2]

# Open up the HDF5 objects
new = h5py.File(new_fname, 'r')
idp = h5py.File(idp_fname, 'r')


def sif_plot(ax1, ax2, ax3, window):

    new_key = f"/physical_retrieval_results/{window}_SIF_absolute"
    idp_key = f"/DOASFluorescence/fluorescence_offset_relative_{window}_idp"

    y = new[new_key][:]
    x = idp[idp_key][:]

    new_chi2 = new[f'/physical_retrieval_results/{window}_retrieved_chi2'][:]
    idp_chi2 = idp[f'/DOASFluorescence/residual_reduced_chi2_fluorescence_{window}_idp'][:]
    mask = np.isfinite(y) & (new_chi2 < 10) & \
        (idp['/DOASFluorescence/fluorescence_qual_flag_idp'][:] == 0)

    y = y[mask] / idp[f'/DOASFluorescence/continuum_level_radiance_{window}_idp'][:][mask]
    x = x[mask]

    lreg = sps.linregress(x, y)

    ex_min = min(x.min(), y.min())
    ex_max = max(x.max(), y.max())

    _x = np.array([ex_min, ex_max])

    ax1.hexbin(x, y,
               mincnt=1, linewidths=0, gridsize=50,
               extent=[ex_min, ex_max, ex_min, ex_max],
               cmap='plasma')
    ax1.set_aspect(1.0)

    ax1.plot(_x, _x, '--', c='grey', linewidth=0.5)
    ax1.plot(_x, _x * lreg[0] + lreg[1], '--', c='red',
             linewidth=0.5, label='Fit')
    ax1.set_ylabel('New SIF (rel)')
    ax1.set_xlabel('IDP SIF (rel)')


    title_str = f"{window}\n"
    title_str += "$y = {:.3f}\cdot x {:+.3f}$\n".format(lreg[0], lreg[1])
    title_str += "$R = {:.3f}$".format(lreg[2])
    ax1.set_title(title_str, fontsize=8)

    ax1.grid(linewidth=0.5, linestyle='dotted')

    ax2.hist(new_chi2[mask], bins=100, range=(0, 2.5), alpha=0.5,
            label="New ({:0.2f})".format(np.median(new_chi2[mask])))
    ax2.hist(idp_chi2[mask], bins=100, range=(0, 2.5), alpha=0.5,
            label="Old ({:0.2f})".format(np.median(idp_chi2[mask])))
    ax2.legend(fontsize=8)
    ax2.set_title("$\chi^2$")

    ax3.hist(new[f'physical_retrieval_results/{window}_num_iterations'][:][mask],
             bins=np.arange(5)-0.5, alpha=0.5, label='New', rwidth=0.35, log=False)
    ax3.hist(idp[f'DOASFluorescence/iterations_fluorescence_{window}_idp'][:][mask],
             bins=np.arange(5)-0.5, alpha=0.5, label='IDP', rwidth=0.35, log=False)
    ax3.legend(fontsize=8)
    ax3.set_title("# Iterations", fontsize=8)

def ratio_plot(ax1, ax2, ax3, gas):

    chi2_weak = new['/physical_retrieval_results/weak_co2_retrieved_chi2'][:]
    chi2_strong = new['/physical_retrieval_results/strong_co2_retrieved_chi2'][:]

    num_iter_weak = new['/physical_retrieval_results/weak_co2_num_iterations'][:]
    num_iter_strong = new['/physical_retrieval_results/strong_co2_num_iterations'][:]


    ratio_new = new[f'physical_retrieval_results/weak_co2_{gas}_scale'][:] /\
        new[f'physical_retrieval_results/strong_co2_{gas}_scale'][:]
    ratio_idp = idp[f'/DOASCloudScreen/{gas.lower()}_ratio_idp'][:]

    mask = (np.isfinite(ratio_new) &
            #(ratio_idp > 0) & (ratio_idp < 3) &
            (idp['/DOASCloudScreen/cloud_flag_idp'][:] > 0) &
            (new[f'physical_retrieval_results/strong_co2_{gas}_scale'][:] > 0.01) &
            (new[f'physical_retrieval_results/strong_co2_{gas}_scale'][:] < 3) &
            (new[f'physical_retrieval_results/weak_co2_{gas}_scale'][:] > 0.01) &
            (new[f'physical_retrieval_results/weak_co2_{gas}_scale'][:] < 3) &
            (chi2_weak < 10) &
            (chi2_strong < 10))

    y = ratio_new[mask]
    x = ratio_idp[mask]

    lreg = sps.linregress(x, y)

    ex_min = 0.98 * min(np.percentile(x,2), np.percentile(y,2))
    ex_max = 1.02 * max(np.percentile(x,98), np.percentile(y,98))

    _x = np.array([ex_min, ex_max])

    ax1.hexbin(x, y, mincnt=1, linewidths=0, gridsize=55,
               extent=[ex_min, ex_max, ex_min, ex_max], cmap='plasma')
    ax1.set_aspect(1.0)

    ax1.plot(_x, _x, '--', c='grey', linewidth=0.5)
    ax1.plot(_x, _x * lreg[0] + lreg[1], '--', c='red',
             linewidth=0.5, label='Fit')

    ax1.grid(linewidth=0.5, linestyle='dotted')

    title_str = "$y = {:.3f}\cdot x {:+.3f}$\n".format(lreg[0], lreg[1])
    title_str += "$R = {:.3f}$".format(lreg[2])
    ax1.set_title(title_str, fontsize=8)

    ax1.set_xlim(ex_min, ex_max)
    ax1.set_ylim(ex_min, ex_max)

    ax1.set_ylabel(f'New {gas} ratio')
    ax1.set_xlabel(f'IDP {gas} ratio')

    chi2_weak = chi2_weak[mask] #[np.isfinite(chi2_weak)]
    chi2_strong = chi2_strong[mask] #[np.isfinite(chi2_strong)]

    chi2_min = min(np.percentile(chi2_weak, 2),
                   np.percentile(chi2_strong, 2)) * 0.95
    chi2_max = max(np.percentile(chi2_weak, 98),
                   np.percentile(chi2_strong, 98)) * 1.05

    ax2.hist(chi2_weak, bins=50,
             range=(chi2_min, chi2_max),
             alpha=0.5, log=False,
             label="WCO2 ({:0.2f})".format(np.median(chi2_weak)))
    ax2.hist(chi2_strong, bins=50,
             range=(chi2_min, chi2_max),
             alpha=0.5, log=False,
             label="SCO2 ({:0.2f})".format(np.median(chi2_strong)))

    ax2.legend(fontsize=8)
    ax2.set_title("$\chi^2$")

    ax3.hist(new[f'physical_retrieval_results/weak_co2_num_iterations'][:][mask],
             bins=np.arange(1,8)-0.5, alpha=0.5, label='WCO2', rwidth=0.35, log=False)
    ax3.hist(new[f'physical_retrieval_results/strong_co2_num_iterations'][:][mask],
             bins=np.arange(1,8)-0.5, alpha=0.5, label='SCO2', rwidth=0.35, log=False)
    ax3.legend(fontsize=8)
    ax3.set_title("# Iterations", fontsize=8)




# Create the figure canvas
fig, axes = plt.subplots(4, 3, figsize=(7, 9)) #, dpi=150)

##########################################################################
# Panel 1 & 2, 757nm SIF

# Do we have SIF in this dataset?
if "/physical_retrieval_results/757nm_SIF_absolute" in new:

    sif_plot(axes[0,0], axes[0,1], axes[0,2], '757nm')
    sif_plot(axes[1,0], axes[1,1], axes[1,2], '771nm')

if (("/physical_retrieval_results/strong_co2_retrieved_chi2" in new) and
    ("/physical_retrieval_results/weak_co2_retrieved_chi2" in new)):

    ratio_plot(axes[2,0], axes[2,1], axes[2,2], 'CO2')
    ratio_plot(axes[3,0], axes[3,1], axes[3,2], 'H2O')


#fig.align_labels()
fig.tight_layout()

plt.savefig('IDP_comparison_new.pdf', bbox_inches='tight')
plt.show()
