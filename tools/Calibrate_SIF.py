import pickle
import numpy as np
import h5py
import sys
from IPython import embed

print("This file: ", sys.argv[1])

l2 = h5py.File(sys.argv[1], 'r+')
calibs = [sys.argv[2], sys.argv[3]]



for i, band in enumerate(['757', '771']):

    bfile = open(calibs[i], 'rb')
    calib_curve = pickle.load(bfile)
    bfile.close()

    sif = l2[f'RetrievalResults/physical/{band}nm/sif_radiance_gbg'][:].copy()
    cont = l2[f'RetrievalResults/physical/{band}nm/continuum_level_radiance_gbg'][:].copy()

    sif_corr = np.zeros_like(sif)
    sif_corr[:] = np.nan

    valid = (cont >= calib_curve.x[0]) & (cont <= calib_curve.x[-1])
    sif_corr[valid] = sif[valid] - calib_curve(cont)[valid]

    if not 'Fluorescence' in l2['HighLevelResults']:
        l2.create_group('/HighLevelResults/Fluorescence')
    else:
        print("Fluorescence group already exists.")

    try:
        if not f'sif_{band}nm' in l2['HighLevelResults/Fluorescence']:
            l2.create_dataset(f"HighLevelResults/Fluorescence/sif_{band}nm",
                              data=sif_corr)
        else:
            l2[f'HighLevelResults/Fluorescence/sif_{band}nm'][:] = sif_corr

        if not f'sif_raw_{band}nm' in l2['HighLevelResults/Fluorescence']:
            l2.create_dataset(f"HighLevelResults/Fluorescence/sif_raw_{band}nm",
                              data=sif)
        else:
            l2[f'HighLevelResults/Fluorescence/sif_raw_{band}nm'][:] = sif
    except:
        print("Could not write SIF results.")
        embed()
