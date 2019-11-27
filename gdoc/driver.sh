#!/bin/bash

control_file_dir="/home/hcronk/geocarb/ditl_1/control/gasbag"

control_files=($(ls ${control_file_dir}))

for control_file in ${control_files[@]}
do
    echo "Processing" $control_file
    OMP_NUM_THREADS=4 /GASBAG/build/GASBAG -c ${control_file_dir}/${control_file}
    #exit
done # control file loop
