# -*- coding: utf-8 -*-
"""ExtractBasisFunctions.py

This script takes a list of L1B files, and extracts all radiance basis
functions for SIF retrievals according to Guanter et al. 2012
(https://doi.org/10.1016/j.rse.2012.02.006).

The measurements over non-vegetation surfaces are identified via a land cover
map; at the moment only the ESA-CCI land cover map is supported.

Here's the workflow of the code:

1.) After checking the command-line arguments for validity, and opening the
    LC data, do a loop over all supplied L1B files.

2.) Match the soundings to approximate LC classes (no exact footprint matching
    is performed to save time), and then select only those where we can safely
    assume no SIF to be present (water, bare, deserts, ice, urban).

3.) Like in Guanter et al., the radiances are slope-normalised, and then some
    additional filtering is done to kick out spectra that deviate from the
    'median' look of the set. Unfortunately, the spike filter in the OCO-2
    L1B does not detect all spikes, and some small ones still remain in the
    typical SIF bands.

4.) Collect all spectra from all files into per-footprint containers, and then
    unleash the SVD computation on them. Singular vectors are then
    'sign-corrected' using the algorithm by Bro et al.

5.) Make some plots of the top SVs and save the data into HDF files. Note that
    the waveforms are saved in 'absolute pixel'-space, so they correspond to
    the radiance arrays in the L1B files. This should make it easier to deal
    with during the retrieval stage.


Runtime is roughly 15 minutes for 10 L1B files, give or take.

"""
from IPython import embed

import argparse
from netCDF4 import Dataset
import numpy as np
import scipy as sp
from scipy import stats, spatial
import sys
import os.path
import logging
import h5py
import palettable
from matplotlib import pyplot as plt

logformat = ("%(asctime)s %(levelname)-8s - "
             "[%(funcName)s:%(lineno)d] %(message)s")
logging.basicConfig(level=logging.DEBUG, format=logformat)
logger = logging.getLogger()


def bro_algorithm(X, u,s,v):

    # Algorithm to flip signs of singular vectors u,v so they all
    # point into the same direction.

    # For details see: https://doi.org/10.1002/cem.1122
    # I think in the first line of Step 3 (Figure 3), it's meant to say
    # signs instead of absolute values - not sure. It seems to do the trick
    # most of the time, though.

    u_new = u.copy()
    v_new = v.copy()

    N_s = len(s)

    for idx_s in range(N_s):
        s_loc = s.copy()
        s_loc[idx_s] = 0.0
        Y = X - np.dot(u * s_loc, v)

        s_left = 0.0
        s_right = 0.0
        for j in range(N_s):
            _tmpdot = np.dot(u[:, idx_s], Y[:, j])
            s_left += np.sign(_tmpdot) * _tmpdot**2
            _tmpdot = np.dot(v[idx_s], Y[j])
            s_right += np.sign(_tmpdot) * _tmpdot**2

        if np.sign(s_left) != np.sign(s_right):
            if np.abs(s_left) < np.abs(s_right):
                s_left = -s_left
            else:
                s_right = -s_right

        u_new[:, idx_s] = np.sign(s_left) * u[:, idx_s]
        v_new[idx_s] = np.sign(s_right) * v[idx_s]

    return u_new, v_new



def filter_radiances_by_similarity(rad_array, percentile=80):

    # Median of all radiances
    median_rad = np.median(rad_array, axis=0)

    # Compute the distance between all measurements and the median
    # radiance in the array.
    dist = np.zeros(rad_array.shape[0])
    for i, spec in enumerate(rad_array):
        dist[i] = sp.spatial.distance.cosine(median_rad, spec)

    # Take out the best X percent (depending on 'percentile') and return the
    # filtered result.
    good = np.where(dist <= np.percentile(dist, percentile))[0]
    return rad_array[good]

def normalize_radiances(rad_array, num_iter=5, max_filter=1.10,
                        fname=''):
    """Here, we normalize an array of radiances (samples, spectral index) by
dividing through a best-fit slope that characterises the continuum

    """

    rad_norm = rad_array.copy()

    # x (wl) coordinate the same for every measurement here
    x = np.arange(rad_array.shape[1])

    # Loop through every measurement
    for m_idx in range(rad_array.shape[0]):
        # Do the slope fit several times
        for fit_iter in range(num_iter):
            # find the continuum level points
            cont_thres = np.percentile(rad_norm[m_idx], 40)
            cont_idx = np.where(rad_norm[m_idx] > cont_thres)[0]
            # perform linear fit
            slope, icept = np.polyfit(x[cont_idx], rad_norm[m_idx, cont_idx], 1)

            rad_norm[m_idx] /= (icept + slope * x)

    # Ok, so the Spike filter for OCO-2 does not seem to catch all of the
    # cosmic rays, so we do some additional filtering here. The basic idea is
    # that the continuum level should be more or less sitting around 1.0 for
    # normalized radiances. If this is not the case, consider the spectrum
    # dropped.
    good_spectra = np.where(rad_norm.max(axis=1) < max_filter)[0]
    perc_good = 100.0 * len(good_spectra) / rad_norm.shape[0]
    logger.info("Keeping {:.2f}% after normalization.".format(perc_good))


    ## temporary - save a nice plot of all considered spectra
    fig = plt.figure()
    plt.plot(rad_norm[good_spectra, :].T, alpha=0.01, c='k', lw=1.0)
    plt.plot(np.max(rad_norm[good_spectra, :], axis=0), lw=1.0, c='blue')
    plt.plot(np.min(rad_norm[good_spectra, :], axis=0), lw=1.0, c='red')
    plt.plot(np.median(rad_norm[good_spectra, :], axis=0), lw=1.0, c='green')
    plt.title(f'N={len(good_spectra)}')
    plt.savefig(fname)
    plt.close()

    return rad_norm[good_spectra, :]

def grab_radiances_oco2(h5, micro_window, bare_soundings,
                        cont_threshold=2e20):
    """This function grabs radinces from the band in question and does some
OCO-2 specific things like checking for spikes etc.

    :param microwindow: tuple/list with ['name', window_start, window_end]
    :returns: dict with numpy array that contains per-footprint radiances
    :rtype: np.array

    """


    disp_coeffs = h5['InstrumentHeader/dispersion_coef_samp'][0,:,:]
    wl_grid = np.zeros(h5['SoundingMeasurements/radiance_o2'].shape[1:])

    mwin_min = micro_window[1]
    mwin_max = micro_window[2]


    # Construct the wavelength grid from the dispersion coefficients
    for fp in range(8):
        num_pix = wl_grid[fp].shape[0]
        wl_grid[fp] = 1000 * np.poly1d(disp_coeffs[fp][::-1])(np.arange(num_pix))

    # Find the pixels which correspond to the chosen microwindow, and store
    # them in a dict for every footprint index (number-1)
    idx_min = dict.fromkeys(range(8))
    idx_max = dict.fromkeys(range(8))

    for fp in range(8):
        # Note! The -1 and +1 are required such that the dispersion limits match
        # up with how it's done in the retrieval algorithm.
        idx_min[fp] = np.searchsorted(wl_grid[fp], mwin_min) - 1
        idx_max[fp] = np.searchsorted(wl_grid[fp], mwin_max) + 1

    # Read all radiances in for faster access:
    all_radiances = h5['SoundingMeasurements/radiance_o2'][:]
    # Read all Spike data:
    all_spikes = h5['SpikeEOF/spike_eof_weighted_residual_o2'][:]

    # Similarly, store all the radiances into arrays for every footprint
    # separately.
    radiance_dict = dict.fromkeys(range(8))
    for fp in range(8):
        # Select for footprint
        idx_fp = bare_soundings[1] == fp
        # And make sure that there are no Spikes in the spectral range
        max_spikes = np.abs(all_spikes[bare_soundings[0][idx_fp], fp,
                                       idx_min[fp]:idx_max[fp]]).max(axis=1)
        # Also keep very dark scenes out
        cont_level = h5['SoundingMeasurements/rad_continuum_o2'][:][bare_soundings[0][idx_fp], fp]

        # And also check for the quality flag
        qual_flag = h5['FootprintGeometry/footprint_o2_qual_flag'][:][bare_soundings[0][idx_fp], fp]
        # Full filter
        full_filt = ((max_spikes == 0) &
                     (cont_level > cont_threshold) &
                     (qual_flag == 0))

        # And select those radiances
        radiance_dict[fp] = all_radiances[bare_soundings[0][idx_fp][full_filt],
                                          fp, idx_min[fp]:idx_max[fp]]

    # Return the data
    return radiance_dict, idx_min, idx_max


def identify_bare_soundings(LC_data, LC_codes,
                            bare_LC_indices, cf_threshold=0.9):
    """From the sampled LC data, we grab the positions of those where the
cumulative fraction of LC indices belonging to 'bare' surfaces (user-defined),
reaches or goes beyond a certain threshold.

    :param LC_data: LC array with LC fractions (shape: frame,footprint,LC)
    :param LC_filepath: Path to the LC file
    :param bare_LC_indices: Which LC classes do we consider 'bare' surfaces?
    :param cf_threshold: What is the cumulative threshold for it to count as
bare surface?
    :returns: index array with positions of eligible soundings
    :rtype: np.array

    """

    # Match the bare_LC_indices with the values in the LC file
    LC_idx = []
    for idx in bare_LC_indices:
        pos = np.where(LC_codes == idx)[0][0]
        LC_idx.append(pos)
    LC_idx = np.array(LC_idx)

    # Sum up the relevant ones to get the cumulative fraction of bare-type
    # LC pixels within the area
    LC_cumulative = LC_data[:,:, LC_idx].sum(axis=-1)

    # And finally determine where the eligible ones are!
    bare_pos = np.where(LC_cumulative > cf_threshold)
    # .. and let the user know how many.
    logger.info(f"We have {len(bare_pos[0])} soundings over bare surfaces.")

    return bare_pos


def grab_LC_data(lon, lat, LC, LC_lon, LC_lat, LC_codes):
    """Grab LC information based on lon, lat of sounding centers. Note that
this is a quick-n-dirty lookup that does not take into account the proper
shape of the footprint; mainly because that would be a bit too costly and slow.
Also for now, only the ESA CCI land cover map is implemented..

    :param lon: Footprint center longitude
    :param lat: Footprint center latitude
    :returns: LC data
    :rtype: np.array

    """


    # Based on LC codes create empty container that will hold our LC
    # coverage data.
    LC_sampled = np.zeros(lon.shape + (len(LC_codes),),  dtype='<f4')

    # First we need to make sure we have valid Lon/Lat pairs, so let's create
    # a new copy and mask those outliers. The outliers will still be computed
    # (sadly), but we can use the masks to flag the results as invalid.
    valid_lon = np.ma.masked_outside(lon, -180, 180)
    valid_lat = np.ma.masked_outside(lat, -90, 90)

    srt_lon = np.searchsorted(LC_lon, valid_lon.flatten())
    # Caution! ESA CCI LC has descending Lat data, so we need to account for
    # this in the matcher.
    srt_lat = len(LC_lat) - 1 - np.searchsorted(LC_lat, valid_lat.flatten(),
                                                sorter=np.argsort(LC_lat))

    srt_lon = np.ma.masked_array(srt_lon, mask=valid_lon.mask.flatten())
    srt_lat = np.ma.masked_array(srt_lat, mask=valid_lat.mask.flatten())

    # From the centers, we go a few indices outwards to get a better idea
    # of the -probable- LC seen by the instrument. Of course this DOES NOT
    # account for LC map distortion, or spacecraft viewing angle.
    # Remember, this is just a solid guess..

    if instrument == 'oco2':
        delta_EW = 6
        delta_NS = 6


    logger.info(f"Processing LC data for N={len(srt_lon)} soundings.")
    for idx in range(len(srt_lon)):

        # Skip masked indices
        if (srt_lon.mask[idx] or srt_lat.mask[idx]):
            continue

        # Check boundaries! We want to grab a section of the LC map like
        # LC[idx_north:idx_south, idx_west:idx_east]
        # For North/South (lats), we don't wrap around, just cut off at poles
        if srt_lat[idx] > delta_NS:
            idx_north = srt_lat[idx] - delta_NS
        else:
            idx_north = 0

        if srt_lat[idx] < (len(LC_lat) - delta_NS):
            idx_south = srt_lat[idx] + delta_NS
        else:
            idx_north = len(LC_lat) - 1

        # For East/West (lons), we just need to make sure that it wraps around
        # for long indices. Indices < 0 wrap around automatically anyway.
        if srt_lon[idx] > len(LC_lon) + delta_EW:
            idx_east = srt_lon[idx] + delta_EW - len(LC_lon)
        else:
            idx_east = srt_lon[idx] + delta_EW

        idx_west = srt_lon[idx] - delta_EW

        # Now grab those LC pixels
        this_LC = LC[idx_north:idx_south, idx_west:idx_east]
        # Get the fractional contribution, first grab unique values and counts
        unique_LC = np.unique(this_LC, return_counts=True)
        percent_LC = unique_LC[1] / len(this_LC.flatten())

        # And finally shove it into the container
        this_idx = np.unravel_index(idx, lon.shape)
        # .. by looping over all unique LC types and putting them into the
        # right place.
        for i, my_LC in enumerate(unique_LC[0]):
            my_LC_pos = np.where(my_LC == LC_codes)[0][0]
            LC_sampled[this_idx][my_LC_pos] = percent_LC[i]

    # Mask out invalid data
    LC_sampled = np.ma.masked_array(LC_sampled, mask=False)
    LC_sampled.mask[valid_lon.mask] = True
    LC_sampled.mask[valid_lat.mask] = True

    logger.info("Finished sampling LC data.")

    # Ideally, the sum of all LC fractions should be 1.0 (give or take)
    N_dodgy = (LC_sampled.sum(axis=2) < 0.999999).sum()
    if (N_dodgy > 0):
        logger.info(f"We seem to have {N_dodgy} potentially "
                    "dodgy LC aggregations.")
    else:
        logger.info("LC aggregation sanity check passed!")

    return LC_sampled



def grab_location(h5, instrument):
    """Grab locations from L1B file handler "h5", for the "instrument". This
function has some instrument-specific branches to deal with the various
field names. (no version-specifics however) Note that we keep the frame/fp
structure here.

    :param h5: HDF5 file handler
    :param instrument: instrument name (str)
    :returns: numpy arrays (lon, lat)
    :rtype: np.array

    """

    if instrument == 'oco2':
        lon = h5['SoundingGeometry/sounding_longitude'][:]
        lat = h5['SoundingGeometry/sounding_latitude'][:]
    else:
        logger.critical(f"{instrument} is not implemented here!")
        sys.exit(1)

    return lon, lat

def check_arguments(args):
    """Checks command-line arguments whether they make sense or not

    :param args: argparse parser object
    :returns: check_status (0 if successful, > 0 if not)
    :rtype: int

    """

    check_status = 0

    # Is the specified instrument supported?
    if args.instrument.lower() not in valid_instruments:
        logger.critical(f"{args.instrument} is not in the list of "
                        "supported instruments: ")
        for instrument in valid_instruments:
            logger.critical(instrument)

        check_status = 1

    for l1b_file in args.l1b.split(','):
        if os.path.isfile(l1b_file) is False:
            logger.critical(f"{l1b_file} is not a valid file!")
            check_status = 1

    return check_status

if __name__ == '__main__':

    # These instruments (or L1B files) are currently supported
    valid_instruments = ['oco2']

    #### Some user settings here ##############################################

    # This value defines what percentile we keep after sorting the
    # slope-corrected spectra by similarity (compared to the median of all
    # considered spectra per L1B file). Lower percentile means lower number
    # that go into the SVD, but also reduces the chance of getting outliers.
    similarity_percentile = 40

    # The max_filter variable is used when filtering the spectra AFTER
    # slope normalization. Ideally, slope-corrected spectra should have their
    # continuum at a value of ~1.0. Any spikes would be then seen with values
    # above that. Any spectrum with values larger than this one is kicked out.
    max_filter = 1.05


    ###########################################################################



    # Read the command-line arguments
    parser = argparse.ArgumentParser(description="GeoCARBSIF Basis Function "
                                     "Extractor")
    parser.add_argument('-i', dest='instrument', required=True,
                        help="Which instrument are we using?")
    parser.add_argument('-l1b', dest='l1b', required=True,
                        help="Path to list of L1B files")
    parser.add_argument('-w', dest='window', required=True,
                        help="Path to window file")
    parser.add_argument('-o', dest='o', required=True,
                        help="Path to output file (will be overwritten)")
    parser.add_argument('-lc', dest='lc', required=True,
                        help="Path to land cover file")
    args = parser.parse_args()

    # First, check if the command-line arguments are sensible
    if check_arguments(args) == 0:
        logger.info("Command-line parameters accepted!")
    else:
        sys.exit(1)

    # Use lowercase for instrument string, makes it a bit more hassle-free..
    instrument = args.instrument.lower()

    # Read the textfile containing the L1B paths
    l1b_list = [x.rstrip() for x in open(args.l1b, 'r').readlines()]

    # Reat the textfile containing the window information
    win_raw = open(args.window, 'r').readlines()
    if len(win_raw) != 1:
        logger.critical("Sorry, cannot have more than one line in window file.")
        sys.exit(1)

    win_split = win_raw[0].rstrip().split(',')
    if len(win_split) != 3:
        logger.critical("Sorry, must have exactly three entries in window file.")
        sys.exit(1)

    window_name = win_split[0]
    window_min = float(win_split[1])
    window_max = float(win_split[2])

    micro_window = (window_name, window_min, window_max)
    ###########################################################################
    ######### LC file section #################################################

    # Read the LC data
    nc = Dataset(args.lc, 'r')
    logger.info(f"Succesfully opened LC data at {args.lc}")

    # And read all necessary data into memory (lon, lat, LC classes).
    logger.info(f"Reading LC data into memory..")

    LC_lon = nc.variables['lon'][:]
    LC_lat = nc.variables['lat'][:]
    #LC = nc.variables['lccs_class']
    # Or do this if you want to read the whole dataset into memory at once
    LC = nc.variables['lccs_class'][:]
    # What are the possible LC codes?
    LC_codes = nc.variables['lccs_class'].flag_values.astype('uint8')

    logger.info(f"Reading done.")

    ###########################################################################
    # Loop through all L1B files supplied by the user and perform the
    # radiance extraction.

    radiance_list = []

    for l1b_file in l1b_list:

        l1b_filename = l1b_file.split('/')[-1]

        # See if we can actually open the L1B file
        h5 = h5py.File(l1b_file, 'r')
        # h5py will throw an OSError if this does not work, so no need for any
        # try/catch here.
        logger.info(f"Opened {l1b_file}")

        # Grab position (lon, lat) from each sounding in the file
        lon, lat = grab_location(h5, instrument)

        # Quick summary of number of soundings
        if instrument == 'oco2':
            # Grab numbers - oco2 is shaped (Frame, Footprint)
            N_fp = lon.shape[1]
            N_frame = lon.shape[0]
            N = N_fp * N_frame

            # And print out for convenience
            if N_fp != 8:
                logger.critical(f"For OCO-2, we expect 8 footprints, but you "
                                f"gave me {N_fp}!")
                sys.exit()

            logger.info(f"We have {N_frame} frames at {N_fp} footprints for a "
                        f"total of {N} soundings.")

        else:
            logger.critical(f"{instrument} is not implemented here!")
            sys.exit(1)

        # Now onto grabbing the LC information
        LC_data = grab_LC_data(lon, lat, LC, LC_lon, LC_lat, LC_codes)

        # And masking those soundings, which are eligible for
        # vegetation-free basis function extraction.
        bare_soundings = identify_bare_soundings(LC_data, LC_codes,
                                                 [190,  # Urban
                                                  200,  # Bare
                                                  201,  # Other bare
                                                  202,  # Also bare
                                                  210,  # Water bodies
                                                  220], # Permanent snow/ice
                                                 cf_threshold=0.95)

        # Grab the radiances that we need, along with pixel indices
        radiances, idx_min, idx_max = \
            grab_radiances_oco2(h5, micro_window, bare_soundings,
                                cont_threshold=1e20)
        # And normalize them w.r.t. the continuum-level slope
        norm_radiances = dict.fromkeys(radiances.keys())
        for fp in range(8):
            # But we also want to filter them afterwards
            temp = normalize_radiances(radiances[fp], num_iter=3,
                                       max_filter=max_filter,
                                       fname=f'{l1b_filename}_FP{fp+1}_allspectra.png')
            norm_radiances[fp] = filter_radiances_by_similarity(
                temp, percentile=similarity_percentile)

            logger.info(f"Keeping {norm_radiances[fp].shape[0]} out of {temp.shape[0]}")
        radiance_list.append(norm_radiances)

    # Assuming the radiance_list dictionaries have a common number of pixels
    # for each footprint, we can stack them into one big array per footprint.
    all_radiances = dict.fromkeys(radiance_list[0])

    for key in all_radiances.keys():
        temp_list = []
        for i in range(len(radiance_list)):
            temp_list.append(radiance_list[i][key])
        # If vstack fails, then it's likely because the dimensions of the
        # sub-arrays don't match, which means the dispersion between orbits
        # has changed. No easy solution here..
        all_radiances[key] = np.vstack(temp_list)


    # Now that we have all radiances, perform the SVD for every FP
    logger.info("Performing SVD..")
    v_vectors = dict.fromkeys(all_radiances)
    s_values = dict.fromkeys(all_radiances)
    u_vectors = dict.fromkeys(all_radiances)

    for key in v_vectors.keys():
        logger.info("Matrix going into SVD: {}".format(shape(all_radiances[key])))
        u, s, v = sp.linalg.svd(all_radiances[key], full_matrices=False)

        # Sign ambiguity! If v is a signular vector, so is -v. We use the
        # algorithm from Bro et al. 2007 to resolve the ambiguity.
        u_new, v_new = bro_algorithm(all_radiances[key], u, s, v)

        v_vectors[key] = v_new
        s_values[key] = s
        u_vectors[key] = u_new
    logger.info("SVD done!")

    # If requested, get some diagnostic plots done so we know a bit about how
    # the result turned out
    diagnostic = True
    if diagnostic:
        logger.info("Producing diagnostic plots")

        # Plot the first n singular vectors using a nice colour scheme, and
        # by each FP.

        n = 8
        colors = palettable.tableau.Tableau_10.mpl_colors
        for fp in v_vectors:
            fig, axarr = plt.subplots(n, 1, figsize=(5, 8), dpi=200)
            for i in range(n):

                ax = axarr[i]
                ax.plot(v_vectors[fp][i], c=colors[i])

                s_var = s_values[fp][i] / s_values[fp].sum() * 100.0
                ax.text(0.98, 0.02,
                        "SV{:d}, Var% = {:0.3f}".format(i+1, s_var),
                        ha='right', va='bottom', fontsize=10,
                        transform=ax.transAxes)

                if ax != (n-1):
                    ax.set_xticklabels([])

            axarr[0].set_title(f"Footprint {fp+1}")
            fig.tight_layout()
            plt.savefig(f"{args.o}_FP{fp+1}_SVs.png", bbox_inches='tight')
            plt.close()


    # And finally stick them into an HDF file where we can use them for
    # retrievals. Old file will be overwritten. Control here how many
    # SVs you actually want in the file.
    N_sv_output = 20
    fname = f'{args.o}_basisfunctions.h5'
    logger.info(f"Writing out the top {N_sv_output} SVDs into the file: "
                f"{fname}")

    with h5py.File(fname, 'w') as h5out:
        for key in v_vectors:
            for i_sv in range(N_sv_output):

                if instrument == 'oco2':
                    outarr = np.zeros(1016)
                    outarr[:] = np.nan

                outarr[idx_min[key]:idx_max[key]] = v_vectors[key][i_sv]
                h5out.create_dataset(f'BasisFunction_SV{i_sv+1}_FP{key+1}',
                                     data=outarr)

    # Close the LC file
    nc.close()
    logger.info("All done. Have a good day.")
