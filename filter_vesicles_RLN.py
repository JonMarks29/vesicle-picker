from vesicle_picker import (
    postprocess,
    helpers,
    external_import,
    external_export
)
from tqdm import tqdm
import argparse
import os

# Read in the job parameters file
parser = argparse.ArgumentParser(description='Part 2/3 for using vesicle picker in relion')
parser.add_argument(
    'parameters',
    type=str,
    help='Path to .ini file containing the parameters for vesicle picking.'
)


args = parser.parse_args()
parameters_filepath = args.parameters
parameters = helpers.read_config(parameters_filepath)

# set IO from parameter file
basedir = parameters.get('io', 'basedir')
motion = parameters.get('io', 'motioncorr')
number_process = parameters.getint('io', 'number_process')

output = os.path.join(parameters.get('io', 'output'), "filter")

if number_process == -1:
    number_process = 10000000


os.makedirs(output, exist_ok=True)

# get list of mics with masks to process
mics = []
for i in os.listdir(os.path.join(basedir,parameters.get('io', 'output'),"masks")):
    mics.append(i)

# restrict processing to specified number of mics
if len(mics) > number_process:
    mics = mics[:number_process]



# Loop over all micrographs in the job directory
for micrograph in tqdm(mics):
#    if micrograph[:-4]+"_vesicles.pkl" in os.listdir(parameters.get('output', 'directory')):
#        continue
#    if micrograph[:-4]+"_vesicles.pkl" not in os.listdir(parameters.get('input', 'directory')):
#        continue

    # Construct the filename of the file to import
    masks_filename = (os.path.join(basedir,parameters.get('io', 'output'),"masks",micrograph)
    )
    
    # Read in the masks from that UID
    masks = external_import.import_masks_from_disk(masks_filename)

    # Filter these masks based on min and max values recorded in job parameters
    filtered_masks = postprocess.apply_filters(masks, parameters_filepath)

    # Don't write out if no masks pass filtering
    if len(filtered_masks) == 0:
        continue


    # Save the vesicles to the output directory
    # on a micrograph-by-micrograph basis
    external_export.export_masks_to_disk(
        filtered_masks,
        os.path.join(output,os.path.splitext(micrograph)[0]+"_filter.pkl"),
        compression='uint64')

