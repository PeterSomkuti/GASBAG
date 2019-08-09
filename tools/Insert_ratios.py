import numpy as np
import h5py
import argparse
import logging
import os.path
import sys

from IPython import embed

logformat = ("%(asctime)s %(levelname)-8s - "
             "[%(funcName)s:%(lineno)d] %(message)s")
logging.basicConfig(level=logging.DEBUG, format=logformat)
logger = logging.getLogger()


def insert_ratios(h5):

    # We need to check if the fields are actually present in the file
    if ("xco2_gbg" not in h5['RetrievalResults/physical/weak_co2/']):
        logger.error("Field: 'xco2_gbg' not found in 'RetrievalResults/physical/weak_co2/'")
        sys.exit(1)
    if ("xco2_gbg" not in h5['RetrievalResults/physical/strong_co2/']):
        logger.error("Field: 'xco2_gbg' not found in 'RetrievalResults/physical/strong_co2/'")
        sys.exit(1)

    for gas in ['co2', 'h2o']:

        # Pick out gas concentration from weak band
        gas_weak = h5[f'RetrievalResults/physical/weak_co2/{gas}_scale_0.00_1.00_gbg'][:]
        gas_weak_ucert = h5[f'RetrievalResults/physical/weak_co2/{gas}_scale_0.00_1.00_uncertainty_gbg'][:]

        # And from strong band
        gas_strong = h5[f'RetrievalResults/physical/strong_co2/{gas}_scale_0.00_1.00_gbg'][:]
        gas_strong_ucert = h5[f'RetrievalResults/physical/strong_co2/{gas}_scale_0.00_1.00_uncertainty_gbg'][:]

        # .. and calculate both ratio and uncertainty on the ratio
        gas_ratio = gas_weak / gas_strong
        gas_ratio_ucert = gas_ratio * np.sqrt((gas_weak_ucert / gas_weak)**2 +
                                              (gas_strong_ucert / gas_strong)**2)

        if ('HighLevelResults' not in h5):
            h5.create_group('HighLevelResults')

        if ('CloudScreen' not in h5['HighLevelResults']):
            h5['HighLevelResults'].create_group('CloudScreen')

        if f'HighLevelResults/CloudScreen/{gas}_ratio' not in h5:
            h5.create_dataset(f'HighLevelResults/CloudScreen/{gas}_ratio', data=gas_ratio)
            logger.info(f"Written out: HighLevelResults/CloudScreen/{gas}_ratio")
            h5.create_dataset(f'HighLevelResults/CloudScreen/{gas}_ratio_error', data=gas_ratio_ucert)
            logger.info(f"Written out: HighLevelResults/CloudScreen/{gas}_ratio_error")
        else:
            logger.info(f'HighLevelResults/CloudScreen/{gas}_ratio already exists. Skipping.')



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Calculate ratio variables")
    parser.add_argument('-i', dest='infile', required=True,
                        help="GASBAG retrieval ouptut file")
    args = parser.parse_args()


    # Check if file exists
    if os.path.isfile(args.infile):
        logger.info(f"File: {args.infile} found.")
    else:
        logger.error(f"File: {args.infile} was not found. Exiting.")
        sys.exit(1)


    # Try opening the file
    try:
        h5 = h5py.File(args.infile, 'r+')
        logger.info(f"Opening {args.infile} was successful.")
    except:
        logger.error(f"Error opening {args.infile}.")
        sys.exit(1)

    # Calculate ratios and stick them into the file
    insert_ratios(h5)
