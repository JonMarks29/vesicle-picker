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
parser = argparse.ArgumentParser(description='Part 3/3 for vesicle picker in relion, generating coordinates')
parser.add_argument(
    'parameters',
    type=str,
    help='Path to .ini file containing the parameters for vesicle picking.'
)
parser.add_argument(
    '--mode',
    type=str,
    default='edge',
    help='Picking mode, surface/edge (default: edge)'
)


args = parser.parse_args()
parameters_filepath = args.parameters
parameters = helpers.read_config(parameters_filepath)

# set IO options from parameter file
basedir = parameters.get('io', 'basedir')
motion = parameters.get('io', 'motioncorr')
number_process = parameters.getint('io', 'number_process')
mode = args.mode

if number_process < 0:
    number_process = 100000000


print("picking in mode:", mode)
print("saving to", os.path.join(parameters.get('io', 'output'), "picks_"+mode))

output = os.path.join(parameters.get('io', 'output'), "picks_"+mode)


os.makedirs(output, exist_ok=True)

# find all mics with filter.pkl files to pass to picking task
mics = []
for i in os.listdir(os.path.join(basedir,parameters.get('io', 'output'),"filter")):
    mics.append(i)

if len(mics) > number_process:
    mics = mics[:number_process]


picks_all = {}
# Loop over all micrographs in the job directory
for micrograph in tqdm(mics[:]):
    # Construct the filename of the file to import
    masks_filename = (
            os.path.join(basedir,parameters.get('io', 'output'),"filter", micrograph)
    )
    

    # Read in the masks from that UID
    masks = external_import.import_masks_from_disk(masks_filename)

    # Apply the erosion or dilation set in the parameters based on surface/edge picking mode
    # and generate picks based on the masks.
    if mode == "surface":
        dilation_radius = int(parameters.get('picking', 'surface_dilation_radius'))
    elif mode == "edge":
        dilation_radius = int(parameters.get('picking', 'edge_dilation_radius'))
    if dilation_radius  < 0:
        masks = postprocess.erode_masks(
        masks,
        erosion=-dilation_radius,
        psize=float(parameters.get('general', 'psize')),
        downsample=int(parameters.get('general', 'downsample'))
    )

    elif dilation_radius > 0:
        masks = postprocess.dilate_masks(
        masks,
        dilation=dilation_radius,
        psize=float(parameters.get('general', 'psize')),
        downsample=int(parameters.get('general', 'downsample'))
    )

    else:
       # Re-find the edges for picking, since they're lost on saving
        masks = [postprocess.find_contour(mask) for mask in masks]



    pick_indices = postprocess.generate_picks(
        masks,
        psize=float(parameters.get('general', 'psize')),
        downsample=int(parameters.get('general', 'downsample')),
        box_size=int(parameters.get('picking', 'box_size')),
        mode=args.mode
    )

    # Generate the picks dataset
    mic_base = micrograph[:-16]
    picks_all[mic_base] = pick_indices 


# write relion style manualpick.star files for import to relion
for mic, picks in picks_all.items():
    os.makedirs(output, exist_ok=True)
    with open(os.path.join(output,mic+"_manualpick.star"), "w") as out:
        out.write('''
data_

loop_
_rlnCoordinateX #1
_rlnCoordinateY #2'''
)

        for ind in range(len(picks[0])):
            out.write("\n")
            out.write(str(picks[0][ind])+"\t"+str(picks[1][ind]))


    


