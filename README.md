# GASBAG

## About and contact

The code is currently maintained by Peter Somkuti at CSU/CIRA (peter.somkuti@colostate.edu), please get in touch regarding bugs, suggestions, etc. 

## Installation instructions

As of now, the code requires a Fortran 2018 compliant compiler, ideally gfortran 8.0 and higher. The CMake script has not been updated to work with the Intel Fortran compiler, however that is on the TODO list. Following software is required for a successful compilation:

* HDF Fortran libraries
* CMake >= 3.5
* BLAS Fortran libraries
* LAPACK Fortran libraries

After pulling the code from the repository, create a build directory (e.g. ``mkdir build``), and from within, execute CMake. Depending on the system, 
it might be necessary to point CMake to the location of various libraries. For HDF5, for example, one would write ``cmake 
-DHDF5_ROOT=/path/to/hdf5/1.10.4 ..``.

The two main build configurations are "debug" and "release", which can be passed to CMake like so: ``cmake -DCMAKE_BUILD_TYPE=Debug ..``. Debug 
enables the usual debug flags "-g -Og" and also displays compiler warnings. The "release" configuration has more aggressive optimization options, as 
well as suppresses any compiler warnings.

## Code strategy and structure

This code is only intended for single-band retrievals (implementing multi-band retrievals would be a bit tough to do), and simple code structure. The 
goal is to have all retrieval settings to be processed on a per-window basis, meaning that for a single execution of the program, N single-band retrievals (called 'window') for the entire L1B file can be processed, with each of the windows having their own exclusive settings. This should make it fairly easy to set up a processing pipeline. The software is driven by one single configuration that uses sections and options via a text-based Apache-style file. Example:

```
[logger]
# This section defines the location of the logfile as well as the log verbosity
logfile = logtest.log
loglevel = 10

[algorithm]
# The algorithm section contains the user-choice of algorithm(s) and related
# options - however this section will most likely move into the windows section.
SIF_algorithm = physical
N_basisfunctions = 8

[instrument]
# The name of the instrument, and any potential instrument-related settings.
name = oco2

[input]
# Inputs are required. These point to the location of the L1B and MET files,
# and any further files that might be required to successfully run the retrievals.
l1b_file = /Users/petersomkuti/Work/OCO-2/oco2_L1bScND_22592a_180930_B8100_181001171532.h5
;/Users/petersomkuti/Work/OCO-2/oco2_L1bScND_22592a_180930_B8100r_181009142216.h5
met_file = /Users/petersomkuti/Work/OCO-2/oco2_L2MetND_22592a_180930_B8100r_181008215027.h5

[solar]
# The solar model is required for a physics-based retrieval algorithm.
solar_file = /Users/petersomkuti/Work/toon_solar/imap_solar.dat
;/Users/petersomkuti/Work/toon_solar/solar_merged_20160127_600_26316_100.out
;/Users/petersomkuti/Work/toon_solar/imap_solar.dat
solar_type = toon

[output]
# The name of the output file that contains the retrieval results.
output_file = output.h5

[window-2]
# Every window defines a single-band retrieval, hence the retrieval settings
# have to be defined for every window separately. This will be extended further
# by the full state vector, and output options etc. The window number merely
# decides the order of execution.

name = 757nm
wl_min = 0.75819
wl_max = 0.75922

albedo_order = 4

dispersion_order = 1
dispersion_perturbation = 1d-6 5d-9 1d-14

;gases = O2
;atmosphere = /Users/petersomkuti/Work/geocarbsif/work/o2_atmosphere.dat

statevector = albedo sif dispersion
dsigma_scale = 5.0

fft_convolution = False

[window-1]
name = 771nm

wl_min = 0.76963
wl_max = 0.77030

albedo_order = 2

dispersion_order = 1
dispersion_perturbation = 1d-6 5d-9 1d-11

gases = O2
atmosphere = /Users/petersomkuti/Work/geocarbsif/work/o2_atmosphere.dat

statevector = albedo sif psurf dispersion
dsigma_scale = 3.0

fft_convolution = True

[gas-1]
# If a retrieval window contains gases, the gasses (cross-referenced by the name) have
# to be defined in a separate section, which defines which gas corresponds to an ABSCO
# file. The number "-1" here has no specific meaning, it's just the code looks for any
# sections named "gas-x", where x=1..99. Maybe we'll add a HITRAN reader too at some point..
name = O2
spectroscopy_type = absco
spectroscopy_file = /Users/petersomkuti/Work/absco/o2_v151005_cia_mlawer_v151005r1_narrow.h5

```

Logging is done through a third-party library (logging by Christopher MacMackin), which lets one assign each log entry an urgency value (debug=10, trivia=20, info=30, warning=40, error=50, fatal=60). Using the config file, one can choose to suppress any log lower than the log-level.

Results are stored in an HDF file where each window-name is used to identify the results in the file.



