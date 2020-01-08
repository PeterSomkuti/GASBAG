#!/bin/bash

#SBATCH --job-name=OCO2_SOLAR
#SBATCH --output=OCO2_testcase_log.out
#SBATCH --error=OCO2_testcase_log.err
#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --cpus-per-task=1
#SBATCH --time=23:59:00
#SBATCH --mem-per-cpu=5000

/home/psomkuti/GASBAG/build/GASBAG -c /home/psomkuti/GASBAG/oco2_template.ini
python /data10/psomkuti/GASBAG/tools/Insert_ratios.py -i oco2_testcase.h5
