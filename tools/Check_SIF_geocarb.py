from matplotlib import rcParams
rcParams['font.family'] = 'monospace'
rcParams['font.serif'] = 'Iosevka Term'

import h5py
import numpy as np
import scipy as sp
from scipy import stats as sps
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import sys
from cartopy import crs as ccrs

gridsize = 1.0

projection = ccrs.Geostationary(central_longitude=-88.0)


def plot_map(ax, lon, lat, data, vmin, vmax, gridsize, cmap, statistic='mean'):

    mask = np.isfinite(data)

    binres = sps.binned_statistic_2d(lon[mask], lat[mask], data[mask],
                                     range=[[-180, 180], [-90, 90]],
                                     bins=[360.0 / gridsize, 180.0/gridsize],
                                     statistic=statistic)

    pcm = ax.pcolormesh(binres[1], binres[2], binres[0].T,
                        cmap=cmap,
                        vmin=vmin, vmax=vmax,
                        transform=ccrs.PlateCarree(),
                        rasterized=True)

    plt.colorbar(pcm, ax=ax, extend='both') 

    cl = ax.coastlines()
    cl.set_rasterized(True)

def plot_hexbin(ax, xdata, ydata):

    ax.hexbin(xdata, ydata, mincnt=1, bins='log',
              gridsize=(80, 30), lw=0)

def plot_hist(ax, data, data_range, errors, bins):

    ax.hist(data, range=(data_range[0], data_range[1]), bins=bins,
            log=False)
    ylims = ax.get_ylim()
    ax.vlines([-np.nanmean(errors), np.nanmean(errors)],
              ymin=ylims[0], ymax=ylims[1],
              linestyle='solid', color='red', lw=1.0)
    ax.vlines([-np.nanstd(data), np.nanstd(data)],
              ymin=ylims[0], ymax=ylims[1],
              linestyle='dashed', color='black', lw=1.0)
    ax.set_ylim(ylims)

####### #######
####### #######
####### #######

plt.ioff()

f = h5py.File(sys.argv[1], 'r')

lon = f['SoundingGeometry/sounding_longitude'][:,0]
lat = f['SoundingGeometry/sounding_latitude'][:,0]

sza = f['SoundingGeometry/sounding_solar_zenith'][:,0]
vza = f['SoundingGeometry/sounding_zenith'][:,0]



if 'physical_retrieval_results' in f:

    land_water = f['SoundingGeometry/sounding_land_water_indicator'][:,0]

    cont_755 = np.nanmax(f['physical_retrieval_results/modelled_radiance_757nm'][:,0], axis=1) / 1e18
    cont_772 = np.nanmax(f['physical_retrieval_results/modelled_radiance_771nm'][:,0], axis=1) / 1e18

    mask_755 = (f['physical_retrieval_results/757nm_num_iterations'][:,0] > 20) & (f['physical_retrieval_results/757nm_retrieved_chi2'][:,0] > 3)
    mask_772 = (f['physical_retrieval_results/771nm_num_iterations'][:,0] > 20) & (f['physical_retrieval_results/771nm_retrieved_chi2'][:,0] > 3)

    fs_755 = f['physical_retrieval_results/757nm_SIF_absolute'][:,0] / 1e18
    fs_755[mask_755] = np.nan

    fs_755_error = f['physical_retrieval_results/757nm_SIF_absolute_uncertainty'][:,0] / 1e18
    fs_755_error[mask_755] = np.nan

    fs_772 = f['physical_retrieval_results/771nm_SIF_absolute'][:,0] / 1e18
    fs_772[mask_772] = np.nan

    fs_772_error = f['physical_retrieval_results/771nm_SIF_absolute_uncertainty'][:,0] / 1e18
    fs_772_error[mask_772] = np.nan

if 'RetrievalResults' in f:

    land_water = np.zeros_like(vza)

    cont_755 = f['RetrievalResults/fluorescence_757nm/Ancillary/continuumLevelRadiance'][:,0] / 1e18
    cont_772 = f['RetrievalResults/fluorescence_771nm/Ancillary/continuumLevelRadiance'][:,0] / 1e18

    mask_755 = (cont_755 > 0)
    mask_772 = (cont_772 > 0)

    '''
    fs_755 = f['RetrievalResults/fluorescence_757nm/RetrievedStateVector/state_vector'][:,0,4]
    fs_772 = f['RetrievalResults/fluorescence_771nm/RetrievedStateVector/state_vector'][:,0,5]
    fs_755 = fs_755 * cont_755 / ( 1 + fs_755)
    fs_772 = fs_772 * cont_772 / ( 1 + fs_772)
    '''

    fs_755 = f['HighLevelResults/Fluorescence/fluorescence_radiance_757nm'][:,0] / 1e18
    fs_755_error = f['HighLevelResults/Fluorescence/fluorescence_radiance_757nm_1sigma'][:,0] / 1e18
    fs_772 = f['HighLevelResults/Fluorescence/fluorescence_radiance_771nm'][:,0] / 1e18
    fs_772_error = f['HighLevelResults/Fluorescence/fluorescence_radiance_771nm_1sigma'][:,0] / 1e18

fs_rel_755 = fs_755 / ( cont_755 - fs_755 )
fs_rel_755[mask_755] = np.nan
fs_rel_755_error = fs_755_error / ( cont_755 - fs_755_error )



fs_rel_772 = fs_772 / ( cont_772 - fs_772 )
fs_rel_772[mask_772] = np.nan
fs_rel_772_error = fs_772_error / ( cont_772 - fs_772_error )

# Check overall statistic, bias and spread

fig = plt.figure(figsize=(20, 12), dpi=150)
gs = gridspec.GridSpec(4, 5, width_ratios=[2, 0.5, 0.5, 1, 1])



ax = plt.subplot(gs[0, 0], projection=projection)

ax.set_title("SIF 755nm Absolute\n[$10^{18}$ ph/s/m2/sr/nm]")
plot_map(ax, lon, lat, fs_755, -3, 3, gridsize, 'BrBG')

ax = plt.subplot(gs[0, 1])
plot_hist(ax, fs_755[land_water == 0], (-3, 8),
          fs_755_error[land_water == 0], 100)
ax.set_title("Land only")

ax = plt.subplot(gs[0, 2])
plot_hist(ax, fs_755[land_water == 1], (-4, 4),
          fs_755_error[land_water == 1], 100)
ax.set_title("Water only")

ax = plt.subplot(gs[0, 3])
plot_hexbin(ax, cont_755[land_water == 0], fs_755[land_water == 0])
ax.set_xlabel("Continuum\n[$10^{18}$ ph/s/m2/sr/nm]")
ax.set_ylabel("SIF absolute\n[$10^{18}$ ph/s/m2/sr/nm]")
ax.set_title("Land only")

ax = plt.subplot(gs[0, 4])
plot_hexbin(ax, cont_755[land_water == 1], fs_755[land_water == 1])
ax.set_xlabel("Continuum\n[$10^{18}$ ph/s/m2/sr/nm]")
ax.set_ylabel("SIF absolute\n[$10^{18}$ ph/s/m2/sr/nm]")
ax.set_title("Water only")



ax = plt.subplot(gs[1, 0], projection=projection)
ax.set_title("SIF 755nm Relative [%]")
plot_map(ax, lon, lat, fs_rel_755 * 100.0, -2, 2, gridsize, 'BrBG')

ax = plt.subplot(gs[1, 1])
plot_hist(ax, fs_rel_755[land_water == 0] * 100, (-3, 5),
          fs_rel_755_error[land_water == 0] * 100, 100)
ax.set_title("Land only")

ax = plt.subplot(gs[1, 2])
plot_hist(ax, fs_rel_755[land_water == 1] * 100, (-4, 4),
          fs_rel_755_error[land_water == 1] * 100, 100)
ax.set_title("Water only")

ax = plt.subplot(gs[1, 3])
plot_hexbin(ax, cont_755[land_water == 0], 100 * fs_rel_755[land_water == 0])
ax.set_xlabel("Continuum\n[$10^{18}$ ph/s/m2/sr/nm]")
ax.set_ylabel("SIF absolute\n[$10^{18}$ ph/s/m2/sr/nm]")
ax.set_title("Land only")

ax = plt.subplot(gs[1, 4])
plot_hexbin(ax, cont_755[land_water == 1], 100 * fs_rel_755[land_water == 1])
ax.set_xlabel("Continuum\n[$10^{18}$ ph/s/m2/sr/nm]")
ax.set_ylabel("SIF relative [%]")
ax.set_title("Water only")




ax = plt.subplot(gs[2, 0], projection=projection)
ax.set_title("SIF 771nm Absolute\n[$10^{18}$ ph/s/m2/sr/nm]")
plot_map(ax, lon, lat, fs_772, -2, 2, gridsize, 'BrBG')

ax = plt.subplot(gs[2, 1])
plot_hist(ax, fs_772[land_water == 0], (-3, 5),
          fs_772_error[land_water == 0], 100)

ax = plt.subplot(gs[2, 2])
plot_hist(ax, fs_772[land_water == 1], (-3, 3),
          fs_772_error[land_water == 1], 100)

ax = plt.subplot(gs[2, 3])
plot_hexbin(ax, cont_772[land_water == 0], fs_772[land_water == 0])
ax.set_xlabel("Continuum\n[$10^{18}$ ph/s/m2/sr/nm]")
ax.set_ylabel("SIF absolute\n[$10^{18}$ ph/s/m2/sr/nm]")
ax.set_title("Land only")

ax = plt.subplot(gs[2, 4])
plot_hexbin(ax, cont_772[land_water == 1], fs_772[land_water == 1])
ax.set_xlabel("Continuum\n[$10^{18}$ ph/s/m2/sr/nm]")
ax.set_ylabel("SIF absolute\n[$10^{18}$ ph/s/m2/sr/nm]")
ax.set_title("Water only")




ax = plt.subplot(gs[3, 0], projection=projection)
ax.set_title("SIF 771nm Relative [%]")
plot_map(ax, lon, lat, fs_rel_772 * 100.0, -1, 1, gridsize, 'BrBG')

ax = plt.subplot(gs[3, 1])
plot_hist(ax, 100 * fs_rel_772[land_water == 0], (-2, 3),
          100 * fs_rel_772_error[land_water == 0], 100)

ax = plt.subplot(gs[3, 2])
plot_hist(ax, 100 * fs_rel_772[land_water == 1], (-3, 5),
          100 * fs_rel_772_error[land_water == 1], 100)

ax = plt.subplot(gs[3, 3])
plot_hexbin(ax, cont_772[land_water == 0], 100 * fs_rel_772[land_water == 0])
ax.set_xlabel("Continuum [$10^{18}$ ph/s/m2/sr/nm]")
ax.set_ylabel("SIF relative [%]")
ax.set_title("Land only")

ax = plt.subplot(gs[3, 4])
plot_hexbin(ax, cont_772[land_water == 1], 100 * fs_rel_772[land_water == 1])
ax.set_xlabel("Continuum [$10^{18}$ ph/s/m2/sr/nm]")
ax.set_ylabel("SIF relative [%]")
ax.set_title("Water only")


plt.tight_layout()
plt.savefig(f.filename+"_SIF_check.pdf")
plt.close('all')

plt.ion()
