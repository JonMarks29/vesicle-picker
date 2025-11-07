from vesicle_picker import (
    preprocess,
    generate_masks,
    postprocess,
    helpers,
    external_import,
    external_export,
    funcs_mrcio
)

import matplotlib.pyplot as plt
import numpy as np
import os
from tqdm import tqdm
import argparse

# Define an import function
def import_mrc(filename):
    """Use funcs_mrcio to open a specified .mrc file"""
    # Read the .mrc file in binary
    micrograph = open(filename, 'rb')

    # Use funcs_mrcio to extract image array and
    # rescale values to lie between [-1, 1]
    image = funcs_mrcio.irdsec_opened(micrograph, 0)

    # Use funcs_mrcio to extract header info
    header = funcs_mrcio.irdhdr_opened(micrograph)

    # Return the rescaled image and header
    return image, header

# get parameter file 
parser = argparse.ArgumentParser(description='Part 2b/3 for running vesicle picker with RELION data, displaying vesicles')
parser.add_argument(
    'parameters',
    type=str,
    help='Path to .ini file containing the parameters for vesicle picking.'
)
parser.add_argument(
    '--n',
    type=int,
    default=10,
    help='Path to picking model'
)


args = parser.parse_args()
parameters_filepath = args.parameters
parameters = helpers.read_config(parameters_filepath)

# IO options from parameter file
basedir = parameters.get('io', 'basedir')
motion = parameters.get('io', 'motioncorr')
number_process = args.n
output = os.path.join(parameters.get('io', 'output'))

os.makedirs(os.path.join(output, "display"), exist_ok=True)

if number_process < 0:
    number_process = 100000000

# find images and restrict to requested number
# by default this is 10 for display
toprocess = os.listdir(os.path.join(basedir, output, "masks"))
if len(toprocess) > number_process:
    toprocess = toprocess[:number_process]


count = 0
for mask in tqdm(toprocess):
    mic_base = mask[:-9]
    micname = mic_base+".mrc"

    # load and preprocess image
    image_fullres, header = import_mrc(os.path.join(basedir,motion,micname))
    image_fullres = image_fullres.astype('float32')

    preprocessed_micrograph = preprocess.preprocess_micrograph(
        image_fullres,
        downsample=parameters.getint('general', 'downsample'),
        lowpass_mode=parameters.get('preprocessing', 'lowpass_mode'),
        d=parameters.getint('preprocessing', 'd'),
        sigmaColor=parameters.getint('preprocessing', 'sigmaColor'),
        sigmaSpace=parameters.getint('preprocessing', 'sigmaSpace'))
 
    # get masks for display
    masks = external_import.import_masks_from_disk(os.path.join(basedir, output,"masks", mask))
    filter_masks = external_import.import_masks_from_disk(os.path.join(basedir, output, "filter",  mic_base+"_mask_filter.pkl"))


    # Create a figure with three subplots arranged in a row
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Plot the first image in the first subplot
    axes[0].imshow(preprocessed_micrograph, cmap="Greys_r")
    axes[0].set_title('Preprocessed Micrograph')

    # Plot the second image in the second subplot
    axes[1].imshow(helpers.sum_masks(masks, 'segmentation'))
    axes[1].set_title('Detected Vesicles')

    # # Plot the third image in the third subplot
    axes[2].imshow(helpers.sum_masks(filter_masks, 'segmentation'))
    axes[2].set_title('Filtered Vesicles')

    # Remove axis labels and ticks
    for ax in axes:
        ax.axis('off')

    plt.savefig(os.path.join(basedir, output, "display", mic_base+"_display.png"))

    count += 1
    if count == 10:
        break

