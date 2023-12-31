import numpy as np
import h5py


# Clear-scene radiance
#f = h5py.File('/data8/ttaylor/data_ttaylor/sim_archive/OCO3_orbit_simulations/r93_g04/clear_fluor-off_surf-MSCM_noise-off/intensity/2012_01/OCO3_sim_r93_intensity_clear_fluor-off_surf-MSCM_noise-off_2012_01.hdf', 'r')
f = h5py.File('/home/gregm/geocarb_testing/data/radiative_transfer/output-20160321_20x10_065w-no_aerosol-brdf_3-with_noise/geocarb_l1b_rx_intensity_20160321_20x10_065w-no_aerosol-brdf_3-with_noise.hdf', 'r')

N = f['SoundingGeometry/sounding_zenith'].shape[0]
SZA = f['SoundingGeometry/sounding_solar_zenith'][:,0]
VZA = f['SoundingGeometry/sounding_zenith'][:,0]
mu_bar = (np.cos(np.deg2rad(SZA)) + np.cos(np.deg2rad(VZA))) / (np.cos(np.deg2rad(SZA)) * np.cos(np.deg2rad(VZA)))

weak_pixels = np.arange(f['SoundingMeasurements/radiance_weak_co2'].shape[2]) + 1
weak_wavelength = np.poly1d(f['InstrumentHeader/dispersion_coef_samp'][1,0][::-1])(weak_pixels)

weak_wl_in = [1.60188, 1.60221, 1.60253]
weak_wl_out = [1.60169, 1.60234, 1.60273]
weak_delta_tau = np.zeros((N, len(weak_wl_in)))
weak_idx_in = np.zeros(len(weak_wl_in), dtype='int')
weak_idx_out = np.zeros(len(weak_wl_in), dtype='int')

for j in range(weak_delta_tau.shape[1]):
    weak_idx_in[j] = np.argmin(np.abs(weak_wavelength - weak_wl_in[j]))
    weak_idx_out[j] = np.argmin(np.abs(weak_wavelength - weak_wl_out[j]))

for j in range(weak_delta_tau.shape[1]):
    log_intensity = np.log(
        f['SoundingMeasurements/radiance_weak_co2'][:, 0, weak_idx_out[j]] /
        f['SoundingMeasurements/radiance_weak_co2'][:, 0, weak_idx_in[j]]
            )

    weak_delta_tau[:,j] = log_intensity / mu_bar

# Print it all out
for j in range(weak_delta_tau.shape[1]):
    print("-----WEAK CO2 -----------")
    print("Wavelength in and out: ", weak_wl_in[j], weak_wl_out[j])
    print("Mean (std) delta tau: {:.2f} ({:.2f})"
          .format(weak_delta_tau[:,j].mean(), weak_delta_tau[:,j].std()))


strong_pixels = np.arange(f['SoundingMeasurements/radiance_strong_co2'].shape[2]) + 1
strong_wavelength = np.poly1d(f['InstrumentHeader/dispersion_coef_samp'][2,0][::-1])(strong_pixels)

strong_wl_in = [2.0555, 2.0673, 2.0751]
strong_wl_out = [2.0559, 2.0677, 2.0753]
strong_delta_tau = np.zeros((N, len(strong_wl_in)))
strong_idx_in = np.zeros(len(strong_wl_in), dtype='int')
strong_idx_out = np.zeros(len(strong_wl_in), dtype='int')

for j in range(strong_delta_tau.shape[1]):
    strong_idx_in[j] = np.argmin(np.abs(strong_wavelength - strong_wl_in[j]))
    strong_idx_out[j] = np.argmin(np.abs(strong_wavelength - strong_wl_out[j]))


for j in range(strong_delta_tau.shape[1]):
    log_intensity = np.log(
        f['SoundingMeasurements/radiance_strong_co2'][:, 0, strong_idx_out[j]] /
        f['SoundingMeasurements/radiance_strong_co2'][:, 0, strong_idx_in[j]]
            )

    strong_delta_tau[:,j] = log_intensity / mu_bar

# Print it all out
for j in range(strong_delta_tau.shape[1]):
    print("------STRONG CO2--------")
    print("Wavelength in and out: ", strong_wl_in[j], strong_wl_out[j])
    print("Mean (std) delta tau: {:.2f} ({:.2f})"
          .format(np.nanmean(strong_delta_tau[:,j]),
                  np.nanstd(strong_delta_tau[:,j]))
          )

ch4_pixels = np.arange(f['SoundingMeasurements/radiance_ch4'].shape[2]) + 1
ch4_wavelength = np.poly1d(f['InstrumentHeader/dispersion_coef_samp'][3,0][::-1])(ch4_pixels)

ch4_wl_in = [2.32194, 2.345, 2.317]
ch4_wl_out = [2.32233, 2.3442, 2.314]
ch4_delta_tau = np.zeros((N, len(ch4_wl_in)))
ch4_idx_in = np.zeros(len(ch4_wl_in), dtype='int')
ch4_idx_out = np.zeros(len(ch4_wl_in), dtype='int')

for j in range(ch4_delta_tau.shape[1]):
    ch4_idx_in[j] = np.argmin(np.abs(ch4_wavelength - ch4_wl_in[j]))
    ch4_idx_out[j] = np.argmin(np.abs(ch4_wavelength - ch4_wl_out[j]))


for j in range(ch4_delta_tau.shape[1]):
    log_intensity = np.log(
        f['SoundingMeasurements/radiance_ch4'][:, 0, ch4_idx_out[j]] /
        f['SoundingMeasurements/radiance_ch4'][:, 0, ch4_idx_in[j]]
            )

    ch4_delta_tau[:,j] = log_intensity / mu_bar

# Print it all out
for j in range(ch4_delta_tau.shape[1]):
    print("------ CH4 --------")
    print("Wavelength in and out: ", ch4_wl_in[j], ch4_wl_out[j])
    print("Mean (std) delta tau: {:.2f} ({:.2f})"
          .format(np.nanmean(ch4_delta_tau[:,j]),
                  np.nanstd(ch4_delta_tau[:,j]))
          )
