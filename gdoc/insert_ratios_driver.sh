#!/bin/bash

eval "$($HOME/miniconda/bin/conda shell.bash hook)"

gasbag_code_dir="/GASBAG"
gasbag_file_dir="/home/hcronk/geocarb/ditl_1/data/L2GSB"

gasbag_files=($(ls ${gasbag_file_dir}))

for gf in ${gasbag_files[@]}
do
    echo "Processing" $gf
    python ${gasbag_code_dir}/tools/Insert_ratios.py -i ${gasbag_file_dir}/${gf}

done #gasbag file loop

