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

f = h5py.File(sys.argv[1], 'r')


def plot_map_and_hist(lon, lat, data, gridsize=1.0,
                      title='', fname='default.pdf'):


    mask = np.isfinite(data)
    # Grid the data to regular lon/lat bins
    binres = sps.binned_statistic_2d(lon[mask], lat[mask], data[mask],
                                     range=[[-180, 180], [-90, 90]],
                                     bins=[360.0 / gridsize, 180.0/gridsize],
                                     statistic='mean')

    data_lower, data_higher = np.nanpercentile(binres[0], [2, 98])


    fig = plt.figure(figsize=(12,4), dpi=150)
    gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1])

    ax = plt.subplot(gs[0], projection=ccrs.Robinson())
    ax.coastlines()

    pcm = ax.pcolormesh(binres[1], binres[2], binres[0].T,
                        transform=ccrs.PlateCarree(),
                        vmin=data_lower, vmax=data_higher,
                        cmap='plasma', rasterized=True)
    ax.set_title("(${:.1f}^\circ$)".format(gridsize))
    plt.colorbar(pcm, ax=ax)

    ax = plt.subplot(gs[1])
    ax.hist(binres[0].flatten(),
            range=(0.95 * data_lower, 1.05 * data_higher),
            bins=20)

    hist_med = np.nanmedian(binres[0].flatten())
    hist_iqr = np.diff(np.nanpercentile(binres[0].flatten(), [25, 75]))[0]

    ax.vlines([hist_med], ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1],
              color='red', linestyle='dashed', lw=0.5)

    text_str = "Median: {:.5e}\nIQR: {:.5e}".format(hist_med, hist_iqr)
    ax.text(0.99, 0.99, text_str, transform=ax.transAxes,
            ha='right', va='top', fontsize=10)
    plt.suptitle(title)
    plt.savefig(fname, bbox_inches='tight')
    plt.close('all')


if 'physical_retrieval_results' in f:

    mask = (f['physical_retrieval_results/weak_co2_num_iterations'][:,:].flatten() < 5)
    mask = mask & (f['physical_retrieval_results/strong_co2_num_iterations'][:,:].flatten() < 5)

    co2_weak = f['physical_retrieval_results/weak_co2_XCO2'][:,:].flatten()[mask] * 1e6
    co2_strong = f['physical_retrieval_results/strong_co2_XCO2'][:,:].flatten()[mask] * 1e6
    co2_ratio = co2_weak / co2_strong

    h2o_weak = f['physical_retrieval_results/weak_co2_XH2O'][:,:].flatten()[mask]
    h2o_strong = f['physical_retrieval_results/strong_co2_XH2O'][:,:].flatten()[mask]
    h2o_ratio = h2o_weak / h2o_strong

    lon = f['SoundingGeometry/sounding_longitude'][:,:].flatten()[mask]
    lat =  f['SoundingGeometry/sounding_latitude'][:,:].flatten()[mask]

if 'HighLevelResults' in f:
    co2_ratio = f['HighLevelResults/CloudScreen/co2_ratio'][:,:].flatten()
    h2o_ratio = f['HighLevelResults/CloudScreen/h2o_ratio'][:,:].flatten()

    co2_weak = f['HighLevelResults/CO2_bands/co2_vcd_weakBand'][:,:].flatten()
    co2_strong = f['HighLevelResults/CO2_bands/co2_vcd_strongBand'][:,:].flatten()

if 'physical_retrieval_results' in f:
    plot_map_and_hist(lon, lat, h2o_strong, gridsize=0.5,
                      title='Strong H$_2$O Scale',
                      fname=f.filename+'_strong_h2o_scale.pdf')

    plot_map_and_hist(lon, lat, h2o_weak, gridsize=0.5,
                      title='Weak H$_2$O Scale',
                      fname=f.filename+'_weak_h2o_scale.pdf')


    for band in ['weak', 'strong']:
        chi2 = f[f'physical_retrieval_results/{band}_co2_retrieved_chi2'][:,:].flatten()[mask]
        plot_map_and_hist(lon, lat,
                          chi2,
                          gridsize=0.5,
                          title=f'{band} $\chi^2$',
                          fname=f.filename+f'_{band}_chi2.pdf')

plot_map_and_hist(lon, lat, co2_weak, gridsize=0.5,
                  title='Weak CO$_2$ Scale',
                  fname=f.filename+'_weak_co2_scale.pdf')

plot_map_and_hist(lon, lat, co2_strong, gridsize=0.5,
                  title='Strong CO$_2$ Scale',
                  fname=f.filename+'_strong_co2_scale.pdf')

plot_map_and_hist(lon, lat, co2_ratio, gridsize=0.5,
                  title='CO$_2$ Ratio',
                  fname=f.filename+'_co2_ratio.pdf')

plot_map_and_hist(lon, lat, h2o_ratio, gridsize=0.5,
                  title='H$_2$O Ratio',
                  fname=f.filename+'_h2o_ratio.pdf')


