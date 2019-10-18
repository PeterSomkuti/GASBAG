import h5py
from netCDF4 import Dataset, Dimension, date2num
import numpy as np
import pandas as pd
import argparse
import logging
import os.path
import sys
from scipy import constants as spc
import Simulator_tools

from IPython import embed

logformat = ("%(asctime)s %(levelname)-8s - "
             "[%(funcName)s:%(lineno)d] %(message)s")
logging.basicConfig(level=logging.DEBUG, format=logformat)
logger = logging.getLogger()





def Populate_Lite_File(h5, nc, logname):

    logger.info("Starting population of LtSIF file..")

    logger.info("Reading log file")
    l1blog = Simulator_tools.read_log(logname)

    # CHANGE ME AFTER DOING POST-PROCESSING


    mask = ((h5['SoundingGeometry/sounding_longitude'][:,:].flatten() != 0) &
            (h5['SoundingGeometry/sounding_latitude'][:,:].flatten() != 0))
    #N_sounding = np.product(h5['SoundingGeometry/sounding_id'].shape)
    N_sounding = mask.sum()
    nc_sounding_dim = nc.createDimension("sounding_dim", N_sounding)

    ###############################
    ## SOLAR ZENITH ANGLE
    ###############################
    nc_sza = nc.createVariable('solar_zenith_angle', 'f4', ('sounding_dim'))
    nc_sza[:] = h5['SoundingGeometry/sounding_solar_zenith'][:,:].flatten()[mask]
    nc_sza.standard_name = 'solar_zenith_angle'
    nc_sza.long_name = "Solar zenith angle"
    nc_sza.unit = "degrees"
    nc_sza.comment = "Solar zenith angle is the angle between the line of sight to the sun and the local vertical"
    logger.info("Written SZA")

    ###############################
    ## VIEWING/SENSOR ZENITH ANGLE
    ###############################
    nc_vza = nc.createVariable('sensor_zenith_angle', 'f4', ('sounding_dim'))
    nc_vza[:] = h5['SoundingGeometry/sounding_zenith'][:,:].flatten()[mask]
    nc_vza.standard_name = 'sensor_zenith_angle'
    nc_vza.long_name = "Sensor zenith angle"
    nc_vza.unit = "degrees"
    nc_vza.comment = "Sensor zenith angle is the angle between the line of sight to the sensor and the local vertical"
    logger.info("Written VZA")

    ###############################
    ## VIEWING/SENSOR AZIMUTH ANGLE
    ###############################
    nc_vaa = nc.createVariable('sensor_azimuth_angle', 'f4', ('sounding_dim'))
    nc_vaa[:] = h5['SoundingGeometry/sounding_azimuth'][:,:].flatten()[mask]
    nc_vaa.standard_name = 'sensor_azimuth_angle'
    nc_vaa.long_name = "Sensor azimuth angle"
    nc_vaa.unit = "degrees"
    nc_vaa.comment = "Azimuth angle between line of sight and local north"
    logger.info("Written VAA")

    ###############################
    ## SOLAR AZIMUTH ANGLE
    ###############################
    nc_saa = nc.createVariable('solar_azimuth_angle', 'f4', ('sounding_dim'))
    nc_saa[:] = h5['SoundingGeometry/sounding_solar_azimuth'][:,:].flatten()[mask]
    nc_saa.standard_name = 'solar_azimuth_angle'
    nc_saa.long_name = "Solar azimuth angle"
    nc_saa.unit = "degrees"
    nc_saa.comment = "Azimuth angle between the solar direction as defined by the sounding location, and the sounding local north"
    logger.info("Written SAA")

    ###############################
    ## TIME/DATE
    ###############################
    date_units = "seconds since 1993-1-1 0:0:0"
    dates = pd.to_datetime(h5['SoundingGeometry/sounding_time_string'][:,:].flatten().astype('str'))
    dates_num = date2num(np.array(dates.tz_localize(None), dtype="object"),
                         units=date_units)

    nc_time = nc.createVariable('time', 'f8', ('sounding_dim'))
    nc_time[:] = dates_num[mask]
    nc_time.unit = date_units
    nc_time.standard_time = "time"
    nc_time.long_name = "time"
    nc_time.comment = "Timestamp (seconds since 1 January 1993)"
    logger.info("Written dates")

    ###############################
    ## LONGITUDE
    ###############################
    nc_lon = nc.createVariable('longitude', 'f4', ('sounding_dim'))
    nc_lon[:] = h5['SoundingGeometry/sounding_longitude'][:,:].flatten()[mask]
    nc_lon.standard_name = 'longitude'
    nc_lon.long_name = "Longitude"
    nc_lon.unit = "degrees_east"
    nc_lon.comment = "Center longitude of the measurement"
    logger.info("Written longitude")

    ###############################
    ## LATITUDE
    ###############################
    nc_lat = nc.createVariable('latitude', 'f4', ('sounding_dim'))
    nc_lat[:] = h5['SoundingGeometry/sounding_latitude'][:,:].flatten()[mask]
    nc_lat.standard_name = 'latitude'
    nc_lat.long_name = "Latitude"
    nc_lat.unit = "degrees_north"
    nc_lat.comment = "Center latitude of the measurement"
    logger.info("Written latitude")

    ###############################
    ## 757 SIF PHYSICAL
    ###############################
    nc_sif757 = nc.createVariable('SIF_757nm_physical', 'f4', ('sounding_dim'))
    nc_sif757[:] = h5['physical_retrieval_results/757nm/SIF_absolute'][:,:].flatten()[mask]
    nc_sif757[:] /= (1.0 / ((spc.h * spc.c / (0.757 * 1e-6))))
    nc_sif757.standard_name = 'toa_outgoing_radiance_per_unit_wavelength_due_to_solar_induced_fluorescence'
    nc_sif757.long_name = "Solar Induced Fluorescence at 757nm"
    nc_sif757.unit = "W/m^2/sr/micron"
    nc_sif757.comment = "Solar induced chlorophyll fluorescence at 757nm"
    logger.info("Written SIF 757 physical")

    #################################
    ## 757 SIF UNCERTAINTIES PHYSICAL
    #################################
    nc_sif757_uc = nc.createVariable('SIF_757nm_physical_uncert', 'f4', ('sounding_dim'))
    nc_sif757_uc[:] = h5['physical_retrieval_results/757nm/SIF_absolute_uncertainty'][:,:].flatten()[mask]
    nc_sif757_uc[:] /= (1.0 / ((spc.h * spc.c / (0.757 * 1e-6))))
    nc_sif757_uc.standard_name = 'sif_757nm_uncert'
    nc_sif757_uc.long_name = "1-sigma uncertainty in retrieved SIF"
    nc_sif757_uc.unit = "W/m^2/sr/micron"
    nc_sif757_uc.comment = "1-sigma statistical uncertainty in solar induced chlorophyll fluorescence at 757nm"
    logger.info("Written 757 SIF physical uncertainties")

    ###############################
    ## 757 SIF DATA-DRIVEN
    ###############################
    nc_sif757 = nc.createVariable('SIF_757nm_datadriven', 'f4', ('sounding_dim'))
    nc_sif757[:] = h5['linear_fluorescence_results/757nm/SIF_absolute'][:,:].flatten()[mask]
    nc_sif757[:] /= (1.0 / ((spc.h * spc.c / (0.757 * 1e-6))))
    nc_sif757.standard_name = 'toa_outgoing_radiance_per_unit_wavelength_due_to_solar_induced_fluorescence'
    nc_sif757.long_name = "Solar Induced Fluorescence at 757nm"
    nc_sif757.unit = "W/m^2/sr/micron"
    nc_sif757.comment = "Solar induced chlorophyll fluorescence at 757nm"
    logger.info("Written 757 data-driven SIF ")

    ####################################
    ## 757 SIF UNCERTAINTIES DATA-DRIVEN
    ####################################
    nc_sif757_uc = nc.createVariable('SIF_757nm_datadriven_uncert', 'f4', ('sounding_dim'))
    nc_sif757_uc[:] = h5['linear_fluorescence_results/757nm/SIF_absolute_uncertainty'][:,:].flatten()[mask]
    nc_sif757_uc[:] /= (1.0 / ((spc.h * spc.c / (0.757 * 1e-6))))
    nc_sif757_uc.standard_name = 'sif_757nm_uncert'
    nc_sif757_uc.long_name = "1-sigma uncertainty in retrieved SIF"
    nc_sif757_uc.unit = "W/m^2/sr/micron"
    nc_sif757_uc.comment = "1-sigma statistical uncertainty in solar induced chlorophyll fluorescence at 757nm"
    logger.info("Written 757 data-driven SIF uncertainties")

    embed()

    ####################################
    ## CLOUD OPTICAL DEPTH
    ####################################
    nc_cloud = nc.createVariable('cloud_od', 'f4', ('sounding_dim'))
    nc_cloud[:] = l1blog['Tau_water_1'] + l1blog['Tau_ice_1']
    nc_cloud.standard_name = 'atmosphere_optical_thickness_due_to_cloud'
    nc_cloud.long_name = "Cloud optical depth (liquid water + ice)"
    nc_cloud.unit = "1"
    nc_cloud.comment = "True cloud optical depth"

    ####################################
    ## LAND FRACTION
    ####################################
    nc_landfrac = nc.createVariable('land_fraction', 'f4', ('sounding_dim'))
    nc_landfrac[:] = h5['SoundingGeometry/sounding_land_fraction'][:,0][mask]
    nc_landfrac.standard_name = 'land_area_fraction'
    nc_landfrac.long_name = "Fraction of footprint that is land"
    nc_landfrac.unit = "percent"
    nc_landfrac.comment = "Fraction of footprint that is land"


if __name__ == '__main__':

    # Read the command-line arguments
    parser = argparse.ArgumentParser(description="GASBAG output to LtSIF")

    parser.add_argument('-i', '--infile',
                        required=True, dest='infile',
                        help="Full path to GASBAG output/retrieval HDF5 file")

    parser.add_argument('-o', '--outfile',
                        required=True, dest='outfile',
                        help="Full path to desired location of LtSIF file")

    parser.add_argument('-l', '--logfile',
                        required=True, dest='logfile',
                        help="Full path to the L1b radiance log file")

    # Parse the command line arguments
    args = parser.parse_args()


    # First, let us check if the retrieval output file exists, and if we can read it
    if os.path.isfile(args.infile):
        logger.info(f"File: {args.infile} found.")
    else:
        logger.error(f"File: {args.infile} was not found. Exiting.")
        sys.exit(1)

    # Let h5py try and open it
    try:
        h5 = h5py.File(args.infile, 'r')
        logger.info(f"h5py successfully opened {args.infile}.")
    except Exception as e:
        logger.error(f"h5py could not open the file {args.infile}.")
        logger.error(f"h5py gave following exception: {e}")
        sys.exit(1)

    # Now check if the directory we want to put the file in actually exists, but
    # this only matters if the outdir does not start with "." or "./"
    if (args.outfile[0] == ".") or (args.outfile[0] == "/"):
        outfile_split = os.path.split(args.outfile)
        if os.path.isdir(outfile_split[0]):
            logger.info(f"Path {outfile_split[0]} exists. Continuing.")
        else:
            logger.error(f"Path {outfile_split[0]} does not exist. Exiting")
            sys.exit(1)

    # So far so good - now let's see if we can create a file in there
    try:
        nc = Dataset(args.outfile, 'w')
        logger.info(f"Successfully created output NC4 file: {args.outfile}")
    except Exception as e:
        logger.error(f"Error creating output NC4 file: {args.outfile}")
        logger.error(f"netCDF4 gave following exception: {e}")
        sys.exit(1)

    # All great! Let's move on and create the lite file
    Populate_Lite_File(h5, nc, args.logfile)
