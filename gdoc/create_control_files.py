from jinja2 import Template
from glob import glob
import os
import sys
import argparse
import re

def write_gasbag_control_file(l1b_file):

    l1b_filename_dict = re.match(L1B_REGEX, os.path.basename(l1b_file)).groupdict()
    
    with open(CONTROL_FILE_TEMPLATE, "r") as tf:
        template = Template(tf.read())

    met_file = "geocarb_L2Met_" + l1b_filename_dict["box"] + "-" + l1b_filename_dict["chunk"] + "_" + l1b_filename_dict["yyyymmdd"] + "hhmm_H00000.h5"
    output_file = "geocarb_" + "L2GSB" + "_" + l1b_filename_dict["box"] + "-" + l1b_filename_dict["chunk"] + "_" + l1b_filename_dict["yyyymmdd"] + "hhmm_H00000.h5"
    subdict = {"h5": "ini", "L2GSB": "L2GSBCtl"}
    control_file = os.path.join(CONTROL_DIR, multiple_replace(subdict, os.path.basename(output_file)))
    log_file = re.sub("h5", "log", os.path.basename(output_file))
    
    render = template.render(log_file=os.path.join(LOG_DIR, log_file), 
                            l1b_file=l1b_file, 
                            met_file=os.path.join(MET_DIR, met_file),
                            output_file=os.path.join(OUTPUT_DIR, output_file))

    with open(control_file, "w") as cf:
        cf.write(render)

    return True


def multiple_replace(dict, text):
  # Create a regular expression  from the dictionary keys
  regex = re.compile("(%s)" % "|".join(map(re.escape, dict.keys())))

  # For each match, look-up corresponding value in dictionary
  return regex.sub(lambda mo: dict[mo.string[mo.start():mo.end()]], text) 


if __name__ == "__main__":
    
    CODE_DIR = os.path.dirname(os.path.realpath(__file__))
    CONTROL_FILE_TEMPLATE = os.path.join(CODE_DIR, "template_control.ini")
    #geocarb_l1b_rx_intensity_20160324_1x1_box3_sa1-with_aerosol-brdf_3_chunk017.h5
    L1B_REGEX = "geocarb_(?P<product>(.*))_rx_intensity_(?P<yyyymmdd>[0-9]{8})_(?P<resolution>(.*))_(?P<box>[boxncsa_0-9]{7,8})-(.*)_(?P<chunk>[chunk0-9]{8}).h5$"
    
    parser = argparse.ArgumentParser(description="GeoCarb GASBAG Control File Generator", prefix_chars="-")
    parser.add_argument("-v", "--verbose", help="Prints some basic information during code execution", action="store_true")
    parser.add_argument("-b", "--l1b_file", help="Full path to L1B input file", nargs="*", default=[])
    parser.add_argument("-d", "--input_dir", help="Path to L1B input files to process", default="/data10/psomkuti/GEOCARB_DAY_IN_LIFE/chunked_files/l1b")
    parser.add_argument("-m", "--met_dir", help="Path to L2Met input files to process", default="/data10/psomkuti/GEOCARB_DAY_IN_LIFE/chunked_files/met")
    parser.add_argument("-c", "--control_dir", help="Path to put GASBAG control files", default="/home/hcronk/geocarb/ditl_1/control/gasbag")
    parser.add_argument("-l", "--log_dir", help="Path to put GASBAG log (goes in control files)", default="/home/hcronk/geocarb/ditl_1/logs/gasbag")
    parser.add_argument("-o", "--output_dir", help="Path to put output GASBAG files (goes in control files)", default="/home/hcronk/geocarb/ditl_1/data/L2GSB")    
    args = parser.parse_args()

    verbose = args.verbose
    l1b_files = args.l1b_file
    input_dir = args.input_dir
    MET_DIR = args.met_dir
    CONTROL_DIR = args.control_dir
    LOG_DIR = args.log_dir
    OUTPUT_DIR = args.output_dir
    
    if not l1b_files:
        if not glob(input_dir):
            print(input_dir + " DNE. Exiting")
            sys.exit()
        l1b_files = glob(os.path.join(input_dir, "*.h5"))      

    if not os.path.exists(CONTROL_DIR):
        print(CONTROL_DIR + " DNE. Creating it.")
        os.makedirs(CONTROL_DIR)        
    
    for l1b_file in l1b_files:
        print(l1b_file)
        if not glob(l1b_file):
            print(l1b_file + " DNE. Moving on")
            continue
        if not re.match(L1B_REGEX, os.path.basename(l1b_file)):
            print(l1b_file + " isn't an L1b file. Expected regex =  " + L1B_REGEX)
            print("Moving on")
            continue
        #try:
        write_gasbag_control_file(l1b_file)
        #except:
        #    print("Problem file!")
        #    continue


#/data10/psomkuti/GEOCARB_DAY_IN_LIFE/chunked_files/l1b/
#geocarb_l1b_rx_intensity_20160324_1x1_box1_sa2-with_aerosol-brdf_3_chunk001.h5
#/home/hcronk/geocarb/ditl_1/data/L2Met
#geocarb_GEOS5_20160324_1x1_box3_sa1-with_aerosol-brdf_3_chunk015.h5
#/home/hcronk/geocarb/ditl_1/control/gasbag
#/home/hcronk/geocarb/ditl_1/data/L2GSB
