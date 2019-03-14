from matplotlib import rcParams
rcParams['font.family'] = 'monospace'
rcParams['font.serif'] = 'Iosevka Term'

import matplotlib as mpl
mpl.use('pdf')

import h5py
import numpy as np
import scipy as sp
from scipy import stats as sps
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import sys
from cartopy import crs as ccrs

from IPython import embed

from sklearn import tree
from sklearn.svm import LinearSVC, SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import Perceptron


# Sampling Density Plots
def plot_sampling_density(lon, lat, projection, gridsize, vmax, filename):

    fig = plt.figure(figsize=(10,4), dpi=300)
    ax = plt.subplot(projection=projection)

    binres = sps.binned_statistic_2d(lon, lat, None,
                                     range=[[-180, 180], [-90, 90]],
                                     bins=[360.0 / gridsize, 180.0/gridsize],
                                     statistic='count')
    plotmat = np.ma.masked_equal(binres[0].T, 0)

    pcm = ax.pcolormesh(binres[1], binres[2], plotmat,
                        transform=ccrs.PlateCarree(),
                        vmin=1, vmax=vmax, cmap='plasma', rasterized=True)
    plt.colorbar(pcm, ax=ax, extend='max')

    ax.coastlines()
    plt.savefig(filename, bbox_inches='tight')
    plt.close('all')

# Function to calculate confusion stats
def calculate_statistics(positive, negative, ref_positive, ref_negative, name=""):

    P = ref_positive.sum()
    N = ref_negative.sum()

    TP = (positive & ref_positive).sum()
    TN = (negative & ref_negative).sum()
    FP = (positive & ref_negative).sum()
    FN = (negative & ref_positive).sum()

    TPR = TP / P
    TNR = TN / N

    FNR = 1 - TPR
    FPR = 1 - TNR

    ACC = (TP + TN) / (P + N)

    FDR = FP / (FP + TP)
    FOR = FN / (FN + TN)

    PPV = TP / (TP + FP)
    NPV = TN / (TN + FN)

    print("\n\n{:s}".format(name))
    print("Total: {:d}".format(len(positive)))
    print("Prediction - Clear: {:d} / Non-clear {:d}".format(positive.sum(), negative.sum()))
    print("Truth      - Clear: {:d} / Non-clear {:d}".format(P, N))
    print("-------------------------------------")
    print("       True positive rate: {:7.1f}% ({:d})".format(TPR * 100, TP))
    print("       True negative rate: {:7.1f}% ({:d})".format(TNR * 100, TN))
    print("      False positive rate: {:7.1f}% ({:d})".format(FPR * 100, FP))
    print("      False negative rate: {:7.1f}% ({:d})".format(FNR * 100, FN))
    print("-------------------------------------")
    print("Positive Predictive Value: {:7.1f}%".format(PPV * 100))
    print("Negative Predictive Value: {:7.1f}%".format(NPV * 100))
    #print("-------------------")
    #print("False discovery rate: {:7.1f}%".format(FDR * 100))
    #print(" False omission rate: {:7.1f}%".format(FOR * 100))
    print("-------------------------------------")
    print("                 Accuracy: {:7.1f}%".format(ACC * 100))
    print("-------------------------------------\n\n")


projection = ccrs.Robinson()

# Read in all the data
new = h5py.File('/data6/psomkuti/OCO3_RUNS/SBAS_OCO3_sim_r93_intensity_cloudy_fluor-on_surf-polBRDF_noise-on_2012_02.h5', 'r')

root = "/data8/ttaylor/data_ttaylor/sim_archive/OCO3_orbit_simulations/r93_g04/cloudy_fluor-on_surf-polBRDF_noise-off/intensity/2012_02"

abp = h5py.File(root+'/oco3_SimABPr154_cloudy_fluor-on_surf-polBRDF_noise-off_2012_02.hdf', 'r')
idp = h5py.File(root+'/oco3_SimIDPb7_cloudy_fluor-on_surf-polBRDF_noise-off_2012_02_absco4-2.nc', 'r')

log = open('/data8/ttaylor/data_ttaylor/sim_archive/OCO3_orbit_simulations/r93_g04/meteorology_and_scene/intensity/2012_02/OCO3_scene_r5_intensity_v05222017_2012_02.log', 'r').readlines()


# Read in the logfile first and turn into usable data
for cnt, line in enumerate(log):
    if "Frame" in line:
        frame_line = line
        frame_line_cnt = cnt
        labels_raw = frame_line.split('  ')
        labels = [x.lstrip().rstrip() for x in labels_raw if x != '']
        labels.pop(0)
        break

dtype = [(x, '<f4') for x in labels[:-2]]

array_list = []
for line in log[frame_line_cnt+4:]:
    array_list.append(line.split()[:-2])

truth_raw = np.array(array_list)

truth = np.empty(truth_raw.shape[0], dtype=dtype)

for cnt, label in enumerate(labels[:-2]):
    truth[label] = truth_raw[:,cnt].astype(truth[label].dtype)

# These are the same for all files (hopefully)
lon = new['SoundingGeometry/sounding_longitude'][:,0]
lat = new['SoundingGeometry/sounding_latitude'][:,0]
sza = new['SoundingGeometry/sounding_solar_zenith'][:,0]
vza = new['SoundingGeometry/sounding_zenith'][:,0]
altitude = new['SoundingGeometry/sounding_altitude'][:,0]
land_water = new['SoundingGeometry/sounding_land_water_indicator'][:,0]


# ABP
dp_abp = abp['ABandRetrieval/surface_pressure_delta_abp'][:,0] / 100
abp_cloudflag = abp['ABandRetrieval/cloud_flag_abp'][:,0]

# Total Scene OD
total_od_1 = np.zeros_like(truth['Frame'])
total_od_1 += truth['Cloud water'] + truth['Cloud ice']



#####################
#####################
CLEAR_THRESHOLD = 1.5
clear_true = (total_od_1 <= CLEAR_THRESHOLD)
#clear_true = (truth['Tau_ice_1'] < 0.1) & (total_od_1 <= CLEAR_THRESHOLD)
nonclear_true = (~clear_true)
#####################
#####################




# Ratios
h2o_ratio = new['physical_retrieval_results/weak_co2/XH2O'][:,0] / new['physical_retrieval_results/strong_co2/XH2O'][:,0]
h2o_ratio_lower = new['physical_retrieval_results/weak_co2/H2O_scale_0.850:1.00'][:,0] / new['physical_retrieval_results/strong_co2/H2O_scale_0.850:1.00'][:,0]
h2o_ratio_mid = new['physical_retrieval_results/weak_co2/H2O_scale_0.00:0.850'][:,0] / new['physical_retrieval_results/strong_co2/H2O_scale_0.00:0.850'][:,0]

co2_ratio = new['physical_retrieval_results/weak_co2/XCO2'][:,0] / new['physical_retrieval_results/strong_co2/XCO2'][:,0]
co2_ratio_mid = new['physical_retrieval_results/weak_co2/CO2_scale_0.300:0.700'][:,0] / new['physical_retrieval_results/strong_co2/CO2_scale_0.300:0.700'][:,0]
co2_ratio_lower = new['physical_retrieval_results/weak_co2/CO2_scale_0.700:1.00'][:,0] / new['physical_retrieval_results/strong_co2/CO2_scale_0.700:1.00'][:,0]

co2_ratio_mid_uncert = co2_ratio_mid * np.sqrt((new['physical_retrieval_results/weak_co2/CO2_scale_0.300:0.700_uncertainty'][:,0] /
                                                new['physical_retrieval_results/weak_co2/CO2_scale_0.300:0.700'][:,0])**2 +
                                               (new['physical_retrieval_results/strong_co2/CO2_scale_0.300:0.700_uncertainty'][:,0] /
                                                new['physical_retrieval_results/strong_co2/CO2_scale_0.300:0.700'][:,0])**2)

co2_ratio_lower_uncert = co2_ratio_lower * np.sqrt((new['physical_retrieval_results/weak_co2/CO2_scale_0.700:1.00_uncertainty'][:,0] /
                                                    new['physical_retrieval_results/weak_co2/CO2_scale_0.700:1.00'][:,0])**2 +
                                                   (new['physical_retrieval_results/strong_co2/CO2_scale_0.700:1.00_uncertainty'][:,0] /
                                                    new['physical_retrieval_results/strong_co2/CO2_scale_0.700:1.00'][:,0])**2)

if ('ZLO' in new['physical_retrieval_results/strong_co2']):
    zlo_strong = new['physical_retrieval_results/strong_co2/ZLO'][:,0]
else:
    zlo_strong = np.zeros_like(co2_ratio)

# Try a simple new threshold-based classification
clear_new = (
    (co2_ratio_lower > 1.00) & (co2_ratio_lower < 1.04) &
    (co2_ratio_mid > 0.8) & (co2_ratio_mid < 1.2) &
    (h2o_ratio_lower > 0.96) & (h2o_ratio_lower < 1.02)
    # (zlo_strong > 0) & (zlo_strong < 1e17)
)

mask_new = (
    np.isfinite(co2_ratio_lower) &
    np.isfinite(co2_ratio_mid) &
    np.isfinite(h2o_ratio) &
    #(land_water == 1) & 
    (new['physical_retrieval_results/strong_co2/converged'][:,0] != -1) &
    (new['physical_retrieval_results/weak_co2/converged'][:,0] != -1)
    )

calculate_statistics(clear_new[mask_new], ~clear_new[mask_new], clear_true[mask_new], nonclear_true[mask_new])

co2_ratio_idp = idp['HighLevelResults/CloudScreen/co2_ratio'][:,0]
h2o_ratio_idp = idp['HighLevelResults/CloudScreen/h2o_ratio'][:,0]

clear_idp = np.zeros_like(land_water).astype(bool)
clear_idp[land_water == 0] = (
    (co2_ratio_idp[land_water == 0] >= 0.9) &
    (co2_ratio_idp[land_water == 0] <= 1.03) &
    (h2o_ratio_idp[land_water == 0] >= 0.88) &
    (h2o_ratio_idp[land_water == 0] <= 1.1)
    )

clear_idp[land_water == 1] = (
    (co2_ratio_idp[land_water == 1] >= 1.0005) &
    (co2_ratio_idp[land_water == 1] <= 1.015) &
    (h2o_ratio_idp[land_water == 1] >= 0.88) &
    (h2o_ratio_idp[land_water == 1] <= 1.03)
    )

mask_idp = (idp['RetrievalResults/strong_co2_band/processing_flag'][:,0] == 0) & \
    (idp['RetrievalResults/weak_co2_band/processing_flag'][:,0] == 0)

# mask_new = mask_new & mask_idp
# mask_idp = mask_new & mask_idp

calculate_statistics(clear_idp[mask_idp], ~clear_idp[mask_idp],
                     clear_true[mask_idp], nonclear_true[mask_idp], "IDP")

plot_sampling_density(lon[mask_idp][clear_idp[mask_idp]],
                      lat[mask_idp][clear_idp[mask_idp]], projection, 2.0, 20, 'IDP_clear.pdf')

calculate_statistics(clear_idp[mask_idp] & (abp_cloudflag[mask_idp] == 0),
                     ~(clear_idp[mask_idp] & (abp_cloudflag[mask_idp] == 0)),
                     clear_true[mask_idp],
                     nonclear_true[mask_idp], "IDP & ABP clear")

calculate_statistics(~((~clear_idp[mask_idp]) | (abp_cloudflag[mask_idp] == 1)),
                     ((~clear_idp[mask_idp]) | (abp_cloudflag[mask_idp] == 1)),
                     clear_true[mask_idp],
                     nonclear_true[mask_idp], "IDP | ABP cloudy")

plot_sampling_density(lon[mask_idp][(clear_idp & (abp_cloudflag == 0))[mask_idp]],
                      lat[mask_idp][(clear_idp & (abp_cloudflag == 0))[mask_idp]],
                      projection, 2.0, 20, 'IDP+ABP_clear.pdf')


# Try a decision tree classifier using all the retrieval results
tree_data = np.vstack([
    co2_ratio_mid,
    #co2_ratio_mid_uncert,
    co2_ratio_lower,
    #co2_ratio_lower_uncert,
    #co2_ratio,
    #h2o_ratio,
    h2o_ratio_mid,
    h2o_ratio_lower,
    zlo_strong,
    #new['physical_retrieval_results/strong_co2/num_iterations'][:,0],
    #land_water,
    #new['physical_retrieval_results/strong_co2/converged'][:,0]
    #new['physical_retrieval_results/weak_co2/SIF_absolute'][:,0],
    #new['physical_retrieval_results/strong_co2/SIF_absolute'][:,0],
    #new['physical_retrieval_results/strong_co2/albedo_order_0'][:,0],
    #new['physical_retrieval_results/weak_co2/albedo_order_0'][:,0],
    dp_abp
    #new['physical_retrieval_results/strong_co2/retrieved_chi2'][:,0],
    #new['physical_retrieval_results/weak_co2/retrieved_chi2'][:,0],
    ]).T

tree_target = np.zeros(len(clear_true))
tree_target[clear_true] = 1
tree_target[nonclear_true] = -1

clear_weight = 1 * clear_true.sum() / len(clear_true)
nonclear_weight = 1.0 - clear_weight

tree_clf = tree.DecisionTreeClassifier(max_depth=3,
                                       #min_samples_split=int(0.05 * sum(mask_new)),
                                       #min_samples_leaf=0.05,
                                       criterion='entropy',
                                       class_weight={-1:nonclear_weight, 1:clear_weight})

tree_fit = tree_clf.fit(tree_data[mask_new], tree_target[mask_new])

try:
    dot_data = tree.export_graphviz(tree_fit, out_file='graph.dot',
                                    feature_names=[
                                        'co2_ratio_mid',
                                        #'co2_ratio_mid_uncert',
                                        'co2_ratio_lower',
                                        #'co2_ratio_lower_uncert',
                                        #'co2_ratio',
                                        #'h2o_ratio',
                                        'h2o_ratio_mid',
                                        'h2o_ratio_lower',
                                        'zlo_strong',
                                        #'num_iterations_sco2',
                                        #'land_water',
                                        #'albedo_weak',
                                        #'albedo_strong',
                                        'dpsurf'
                                        #'weak_chi2',
                                        #'strong_chi2'
                                ],
                                    class_names=['nonclear', 'clear'],
                                    filled=True, rounded=True,
                                    leaves_parallel=True, rotate=True)
except:
    print("Failed exporting tree")


tree_predict = tree_fit.predict(tree_data[mask_new])
tree_clear = tree_predict == 1
tree_nonclear = tree_predict == -1

calculate_statistics(tree_clear, tree_nonclear,
                     clear_true[mask_new], nonclear_true[mask_new],
                     "Decision tree")

plot_sampling_density(lon[mask_new][tree_clear], lat[mask_new][tree_clear], projection, 2.0, 20, 'Tree_clear.pdf')

calculate_statistics(tree_clear & (abp_cloudflag[mask_new] == 0),
                     ~(tree_clear & (abp_cloudflag[mask_new] == 0)),
                     clear_true[mask_new], nonclear_true[mask_new],
                     "Decision tree clear & ABP clear")

calculate_statistics(~(tree_nonclear & (abp_cloudflag[mask_new] == 1)),
                     (tree_nonclear & (abp_cloudflag[mask_new] == 1)),
                     clear_true[mask_new], nonclear_true[mask_new],
                     "Decision tree nonclear & ABP cloudy")

calculate_statistics(abp_cloudflag == 0,
                     abp_cloudflag == 1,
                     clear_true[abp_cloudflag != -1],
                     nonclear_true[abp_cloudflag != -1], "ABP")


# Look at OD clear!
fig = plt.figure(figsize=(10, 6))
plt.hist(total_od_1, alpha=0.5, bins=200, range=(0, 7.5),
         label="All", log=True, histtype='step', lw=2.0);
plt.hist(total_od_1[mask_new][tree_clear], alpha=0.5, bins=200,
         range=(0, 7.5), label="Decision Tree clear only", log=True, histtype='step');
plt.hist(total_od_1[mask_new][tree_clear & (abp_cloudflag[mask_new] == 0)], alpha=0.5, bins=200,
         range=(0, 7.5), label="Decision Tree clear + ABP clear", log=True, histtype='step');
plt.hist(total_od_1[clear_idp], alpha=0.5, bins=200,
         range=(0, 7.5), label="IDP clear only", log=True, histtype='step');
plt.hist(total_od_1[abp_cloudflag == 0], alpha=0.5, bins=200,
         range=(0, 7.5), label="ABP clear only", log=True, histtype='step');
plt.ylim(1, 10000)
plt.vlines([0.4], ymin=1, ymax=10000)
plt.xlabel("Total $\\tau$ in Band 1")
plt.legend()
plt.savefig("OD_histogram_clear.pdf", bbox_inches='tight')

# Look at OD cloudy
fig = plt.figure(figsize=(10, 6))
plt.hist(total_od_1, alpha=0.5, bins=200, range=(0, 7.5),
         label="All", log=True, histtype='step', lw=2.0);
plt.hist(total_od_1[mask_new][~tree_clear], alpha=0.5, bins=200,
         range=(0, 7.5), label="Decision Tree nonclear only", log=True, histtype='step');
plt.hist(total_od_1[mask_new][(~tree_clear) | (abp_cloudflag[mask_new] == 1)], alpha=0.5, bins=200,
         range=(0, 7.5), label="Decision Tree nonclear + ABP nonclear", log=True, histtype='step');
plt.hist(total_od_1[~clear_idp], alpha=0.5, bins=200,
         range=(0, 7.5), label="IDP nonclear only", log=True, histtype='step');
plt.hist(total_od_1[abp_cloudflag == 1], alpha=0.5, bins=200,
         range=(0, 7.5), label="ABP nonclear only", log=True, histtype='step');
plt.ylim(1, 10000)
plt.vlines([0.4], ymin=1, ymax=10000)
plt.xlabel("Total $\\tau$ in Band 1")
plt.legend()
plt.savefig("OD_histogram_nonclear.pdf", bbox_inches='tight')

'''
neigh = SVC(gamma=0.001)
nfit = neigh.fit(tree_data[mask_new], tree_target[mask_new])
neigh_predict = nfit.predict(tree_data[mask_new])
neigh_clear = neigh_predict == 1
neigh_nonclear = neigh_predict == -1

calculate_statistics(neigh_clear, neigh_nonclear,
                     clear_true[mask_new], nonclear_true[mask_new],
                     "KNeigh")
'''



for od_type in ['water', 'ice']: #, 'aerosol', 'total']:
    if od_type == 'total':
        this_od = total_od_1
    else:
        this_od = truth[f'Cloud {od_type}']

    fig = plt.figure(figsize=(10, 4))
    lbins = np.logspace(-2, 1, 40)

    abp_clear = abp_cloudflag.copy()
    abp_clear[abp_cloudflag != 0] = 0
    abp_clear[abp_cloudflag == 0] = 1

    binres = sps.binned_statistic(this_od, abp_clear,
                                  statistic='mean', bins=lbins)
    binres2 = sps.binned_statistic(this_od[mask_new],
                                   clear_new[mask_new].astype(int),
                                   statistic='mean', bins=lbins)
    binres3 = sps.binned_statistic(this_od[mask_new],
                                   (clear_new[mask_new] & (abp_cloudflag[mask_new] == 0)).astype(int),
                                   statistic='mean', bins=lbins)
    binres4 = sps.binned_statistic(this_od[mask_idp],
                                   clear_idp[mask_idp].astype(int),
                                   statistic='mean', bins=lbins)
    binres5 = sps.binned_statistic(this_od[mask_new],
                                   tree_clear.astype(int),
                                   statistic='mean', bins=lbins)

    ax1 = plt.gca()

    ax1.semilogx((binres[1][1:] + binres[1][:-1]) * 0.5, binres[0],
                 label='ABP Clear ({:.1f}%)'.format(100 * (abp_cloudflag == 0).sum() / len(abp_cloudflag)),
                 color='blue')
    ax1.semilogx((binres4[1][1:] + binres4[1][:-1]) * 0.5, binres4[0],
                 label='IDP Clear ({:.1f}%)'.format(100 * clear_idp[mask_idp].sum() / len(clear_idp[mask_idp])),
                 color='orange')
    ax1.semilogx((binres2[1][1:] + binres2[1][:-1]) * 0.5, binres2[0],
                 label='Eyeballing Clear ({:.1f}%)'.format(100 * clear_new.sum() / len(clear_new)),
                 color='lightgreen')
    ax1.semilogx((binres5[1][1:] + binres5[1][:-1]) * 0.5, binres5[0],
                 label='Decisition Tree Clear ({:.1f}%)'.format(100 * tree_clear.sum() / len(tree_clear)),
                 color='green')
    ax1.semilogx((binres3[1][1:] + binres3[1][:-1]) * 0.5, binres3[0],
                 label='ABP Clear & Eyeballing Clear ({:.1f}%)'
                 .format(100 * (clear_new[mask_new] & (abp_cloudflag[mask_new] == 0)).sum() / len(tree_clear)),
                 color='black')

    ax2 = plt.twinx()

    ax2.hist(this_od, bins=lbins, alpha=0.85, log=True, zorder=-99, color='grey')

    ax1.set_zorder(ax2.get_zorder()+1)
    ax1.patch.set_visible(False) 
    ax1.legend()
    ax1.set_xlabel("Total $\\tau$")
    ax1.set_ylabel("Fraction Clear")
    if od_type == 'total':
        ax1.vlines([0.4], ymin=0, ymax=1, linestyle='dashed')
    ax2.set_ylabel("Total Scenes")
    ax1.set_title(od_type.capitalize())
    ax2.set_ylim(1, 1e5)

    plt.savefig(f"{od_type}_tommyplot.pdf", bbox_inches='tight')

