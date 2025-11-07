from vesicle_picker import (
    preprocess,
    generate_masks,
    postprocess,
    helpers,
    external_import,
    external_export,
    funcs_mrcio
)
from tqdm import tqdm
import argparse
import os

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



# Read in the job parameters file
parser = argparse.ArgumentParser(description='Part 1/3 for running vesicle picker with RELION data, finding vesicles')
parser.add_argument(
    'parameters',
    type=str,
    help='Path to .ini file containing the parameters for vesicle picking.'
)


args = parser.parse_args()
parameters_filepath = args.parameters
parameters = helpers.read_config(parameters_filepath)

# set IO options from parameter file
basedir = parameters.get('io', 'basedir')
motion = parameters.get('io', 'motioncorr')
number_process = parameters.getint('io', 'number_process')
output = os.path.join(parameters.get('io', 'output'), "masks")

model = parameters.get('segmentation', 'model_weights_path')

if number_process == -1:
    number_process = 10000000



# Initialize the model
model = generate_masks.initialize_model(
    model_weights_path=model,
    model_type=parameters.get('segmentation', 'model_type'),
    device=parameters.get('segmentation', 'device')
)

os.makedirs(output, exist_ok=True)

# check for already processed images
mics = []
done = [] 
for i in os.listdir(os.path.join(basedir,output)):
    done.append(i[:-9])

# find all motioncorr mics
for i in os.listdir(os.path.join(basedir, motion)):
    if "PS" not in i and ".mrc" in i and os.path.splitext(i)[0] not in done:
        mics.append(i)

# only do as many micrographs as requested
if number_process < len(mics):
    mics = mics[:number_process]


# Loop over all micrographs
for micrograph in tqdm(mics):

    # Extract the image
    image_fullres, header = import_mrc(os.path.join(basedir,motion,micrograph))
    image_fullres = image_fullres.astype('float32')

    # Use the preprocess module to get micrograph ready for segmentation.
    # This script uses bilateral filtering,
    # can be adjusted for Gaussian if desired.
    preprocessed_micrograph = preprocess.preprocess_micrograph(
        image_fullres,
        downsample=parameters.getint('general', 'downsample'),
        lowpass_mode=parameters.get('preprocessing', 'lowpass_mode'),
        d=parameters.getint('preprocessing', 'd'),
        sigmaColor=parameters.getint('preprocessing', 'sigmaColor'),
        sigmaSpace=parameters.getint('preprocessing', 'sigmaSpace'))

    # Generate masks with user-optimized parameters
    masks = generate_masks.generate_masks(
        preprocessed_micrograph,
        model,
        psize=parameters.getfloat('general', 'psize'),
        downsample=parameters.getint('general', 'downsample'),
        points_per_side=parameters.getint('segmentation', 'points_per_side'),
        points_per_batch=parameters.getint('segmentation', 'points_per_batch'),
        pred_iou_thresh=parameters.getfloat('segmentation', 'pred_iou_thresh'),
        stability_score_thresh=(
            parameters.getfloat('segmentation', 'stability_score_thresh')
        ),
        crop_n_layers=parameters.getint('segmentation', 'crop_n_layers'),
        crop_n_points_downscale_factor=(
            parameters.getint('segmentation', 'crop_n_points_downscale_factor')
        ),
        crop_nms_thresh=parameters.getfloat('segmentation', 'crop_nms_thresh'),
        min_mask_region_area=(
            parameters.getint('segmentation', 'min_mask_region_area')
        )
    )

    if len(masks) == 0:
        continue

    # Use the postprocess module to compute statistics
    # on the vesicles for downstream filtering
    postprocessed_masks = postprocess.postprocess_masks(
        masks,
        eval(parameters.get('postprocessing', 'functions')),
        preprocessed_micrograph
    )

    # Modify the ellipse fitting key-value pairs to convert to Angstrom
    for mask in postprocessed_masks:
        if 'average_radius' in mask:
            psize = parameters.getfloat('general', 'psize')
            downsample = parameters.getint('general', 'downsample')
            mask['average_radius_A'] = (
                mask['average_radius'] * psize * downsample
            )
            mask['semi_minor_A'] = (
                mask['semi_minor'] * psize * downsample
            )
            mask['semi_major_A'] = (
                mask['semi_major'] * psize * downsample
            )


    # Save the vesicles to the output directory
    # on a micrograph-by-micrograph basis
    external_export.export_masks_to_disk(
        postprocessed_masks,
        os.path.join(basedir,output,os.path.splitext(micrograph)[0]+"_mask.pkl"),
        compression='uint64'
        )
