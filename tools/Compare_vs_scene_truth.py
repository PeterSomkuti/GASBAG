import numpy as np
import h5py
from scipy import interpolate as spi
from scipy import stats as sps
from cartopy import crs as ccrs
from matplotlib import pyplot as plt

prior = np.genfromtxt("ch4_atmosphere_short.dat")
l2 = h5py.File('/home/codell/geocarb/sims/20190409/65Wv2/with_noise-ch4_co_profiles-with_derivs/merged/20160321_20x10_065w-no_aerosol-brdf_3-with_noise-ch4_co_profiles-with_derivs_small_merged.h5', 'r')
scene = h5py.File("/home/gregm/geocarb_testing/data/scene_definition/output-20160321_20x10_065w/geocarb_scene_rx_intensity_20160321_20x10_065w.hdf", 'r')
gbg = h5py.File("ch4_test12.h5", 'r')

land = gbg['SoundingGeometry/sounding_land_fraction'][:,0]
intersect = np.intersect1d(l2['Sound/sounding_id'][:], gbg['SoundingGeometry/sounding_id'][:,0])
srt_scene = np.searchsorted(gbg['SoundingGeometry/sounding_id'][:,0], intersect)
srt_l2 = np.searchsorted(gbg['SoundingGeometry/sounding_id'][:,0], l2['Sound/sounding_id'][:])

xch4_layer = scene['Simulation/Gas/species_density'][:,6,:] / scene['Simulation/Gas/species_density'][:,1,:]
xch4_layer[xch4_layer[:,:] == 0] = np.nan

xch4_truth = scene['Simulation/Gas/species_density'][:,6,:].sum(axis=1) / scene['Simulation/Gas/species_density'][:,1,:].sum(axis=1)
xch4_truth3 = np.zeros_like(xch4_truth) * np.nan
xch4_truth3[srt_l2] = l2['Truth/xch4'][:]

xch4_gbg = gbg['RetrievalResults/physical/ch4/xch4_gbg'][:,0]
xch4_gbg[land == 0] = np.nan

xch4_level = np.zeros((xch4_layer.shape[0], 26))
for i in range(26):

    if (i == 0):
        xch4_level[:,i] = xch4_layer[:,i]
    elif (i == 25):
        xch4_level[:,i] = xch4_layer[:,i-1]
    else:
        xch4_level[:,i] = 0.5 * (xch4_layer[:,i-1] + xch4_layer[:,i])

pressure_levels = scene['Simulation/Thermodynamic/pressure_level'][:]
pressure_levels[:,0] = 1.0
ret_psurf = gbg['RetrievalResults/physical/ch4/surface_pressure_apriori_gbg'][:,0]
xch4_sampled = np.zeros((xch4_layer.shape[0], prior.shape[0]-1))*np.nan
xch4_prior = np.zeros((xch4_layer.shape[0], prior.shape[0]-1))*np.nan

xch4_col_ak = gbg['RetrievalResults/physical/ch4/xch4_column_ak_gbg'][:,:,0].T
xch4_pwgts = gbg['RetrievalResults/physical/ch4/xch4_pressure_weights_gbg'][:,:,0].T

for i in range(xch4_level.shape[0]):
    print(i)
    surf_level = np.searchsorted(prior[1:,0], ret_psurf[i])
    this_p = prior[1:,0][:surf_level + 1].copy()
    this_p[-1] = ret_psurf[i]

    xch4_sampled[i, :surf_level+1] = np.interp(np.log(this_p),
                                               np.log(pressure_levels[i,:26]),
                                               xch4_level[i])

    xch4_prior[i, :surf_level+1] = np.interp(np.log(this_p),
                                             np.log(prior[1:,0]),
                                             prior[1:,4])


xch4_ak_corr = np.nansum((xch4_col_ak - xch4_pwgts) * (xch4_sampled - xch4_prior), axis=1)
xch4_truth2 = np.nansum(xch4_sampled * xch4_pwgts, axis=1)

bins = 1.0
alt = gbg['SoundingGeometry/sounding_altitude'][:,0]
lon = gbg['SoundingGeometry/sounding_longitude'][:,0]
lat = gbg['SoundingGeometry/sounding_latitude'][:,0]
sza = gbg['SoundingGeometry/sounding_solar_zenith'][:,0]


#xch4_gbg -= np.nanmean(xch4_gbg - xch4_truth3/1e9)

xch4_gridded = sps.binned_statistic_2d(
    lon, lat, (xch4_gbg - xch4_ak_corr) * 1e9,
    bins=[360/bins, 180/bins],
    range=[[-180 ,180], [-90, 90]])

dxch4_gridded = sps.binned_statistic_2d(
    lon, lat, (xch4_gbg - xch4_ak_corr - xch4_truth3/1e9) * 1e9,
    bins=[360/bins, 180/bins],
    range=[[-180 ,180], [-90, 90]])


proj = ccrs.Geostationary(central_longitude=-65)
trans = ccrs.PlateCarree()

fig = plt.figure(figsize=(14, 4), dpi=200)

ax = fig.add_subplot(1, 3, 1, projection=proj)
ax.coastlines('50m', lw=0.25)
pcm = ax.pcolormesh(xch4_gridded[1], xch4_gridded[2], xch4_gridded[0].T,
                    cmap='viridis', vmin=1800, vmax=1900, transform=trans)
plt.colorbar(pcm, ax=ax, label='Retrieved XCH$_4$ [ppb]', extend='both')
ax.set_global()

ax = fig.add_subplot(1, 3, 2, projection=proj)
ax.coastlines('50m', lw=0.25)
pcm = ax.pcolormesh(dxch4_gridded[1], dxch4_gridded[2], dxch4_gridded[0].T,
                    cmap='RdYlBu_r', vmin=-25, vmax=25, transform=trans)
plt.colorbar(pcm, ax=ax, label='XCH$_4$ Error [ppb]', extend='both')
ax.set_global()

ax = fig.add_subplot(1, 3, 3)
delta = (xch4_gbg - xch4_ak_corr - xch4_truth3/1e9) * 1e9
ax.hist(delta, bins=100, range=(-50, 50), histtype='step')
ax.text(0.95, 0.95,
        "$\Delta$ = {:.1f} ppb\n$\sigma$ = {:.1f} ppb"
        .format(np.nanmean(delta), np.nanstd(delta)),
        ha='right', va='top', transform=ax.transAxes)
ax.set_xlabel("XCH$_4$ Error [ppb]")
plt.tight_layout()
plt.savefig("XCH4_ak_correct.png", bbox_inches='tight')
plt.close()
