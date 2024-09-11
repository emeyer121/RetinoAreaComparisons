import os
from pprint import pprint
import matplotlib.pyplot as plt

import nibabel as nib
import numpy as np
from nilearn.plotting import plot_roi
from nilearn.image import concat_imgs, math_img, binarize_img
from nilearn.image import resample_img

from nimare.dataset import Dataset
from nimare.decode import discrete
from nimare.utils import get_resource_path
from nimare.extract import download_abstracts, fetch_neuroquery, fetch_neurosynth
from nimare.io import convert_neurosynth_to_dataset

from scipy.sparse import csr_matrix

import logging
import random


def combine_rois(nii_dir):
    # Define the ROI mappings
    roi_to_region = {
        'V1': 'occipital', 'V2': 'occipital', 'V3': 'occipital', 'hV4': 'occipital',
        'VO1': 'ventral', 'VO2': 'ventral', 'PHC1': 'ventral', 'PHC2': 'ventral',
        'V3A': 'dorsal', 'V3B': 'dorsal',
        'LO1': 'lateral', 'LO2': 'lateral', 'TO1': 'lateral', 'TO2': 'lateral',
        'IPS0': 'parietal', 'IPS1': 'parietal', 'IPS2': 'parietal', 'IPS3': 'parietal', 'IPS4': 'parietal', 'IPS5': 'parietal'
    }

    # Initialize the ROI groups
    roi_groups = {
        'occipital': [],
        'ventral': [],
        'dorsal': [],
        'lateral': [],
        'parietal': []
    }

        # Iterate through each file in the directory
    for filename in os.listdir(nii_dir):
        if filename.endswith('.nii') or filename.endswith('.nii.gz'):
            # Load the NIfTI file
            nii_file_path = os.path.join(nii_dir, filename)
            nii_img = nib.load(nii_file_path)
            
            # Extract the ROI name from the filename (assuming the filename contains the ROI name)
            roi_name = os.path.splitext(filename)[0][:-3]
            # Map the ROI to its region
            region = roi_to_region.get(roi_name, None)
            if region:
                roi_groups[region].append(nii_img)

    # Combine the ROIs into regions
    combined_rois = {}
    for region, imgs in roi_groups.items():
        if imgs:
            combined_img = concat_imgs(imgs)
            combined_img = math_img('np.sum(img, axis=-1)', img=combined_img)
            combined_rois[region] = combined_img

    return combined_rois

def get_neurosynth_db():
    out_dir = os.path.abspath("../example_data/")
    os.makedirs(out_dir, exist_ok=True)

    files = fetch_neurosynth(
        data_dir=out_dir,
        version="7",
        overwrite=False,
        source="abstract",
        vocab="terms",
    )
    # Note that the files are saved to a new folder within "out_dir" named "neurosynth".
    neurosynth_db = files[0]

    return neurosynth_db

# Extract region ROIS
nii_dir = os.path.join(os.path.dirname(__file__), 'rois')
combined_rois = combine_rois(nii_dir)

# Create a figure with subplots
regions = list(combined_rois.keys())
# fig, axes = plt.subplots(len(regions), 1, figsize=(10, len(regions) * 5))
count = 0
cmap = ['Blues','Reds','Oranges','Greens','BuGn_r']

for region, data in combined_rois.items():
    print(region)
    # Get the image data array
    data_array = data.get_fdata()
    
    # Create a NIfTI image from the chunked data
    img = nib.Nifti1Image(data_array, affine=data.affine, header=data.header)
    img_bin = binarize_img(img)
    # plt.figure(figsize=(10, 5))
    # plot_roi(img_bin, draw_cross=False,cmap=cmap[count], alpha=1)
    # plt.savefig('./rois/regionROI_'+region+'.png',dpi=1200)
    # plt.savefig('./rois/regionROI_'+region+'.svg',dpi=1200)
    # plt.show()
    count += 1

expansion_factors = [1.3, 2.67, 3.94, 4.25, 7.11]

# Initialize an empty array for the combined data
combined_data = None
count = 0

for region, data in combined_rois.items():
    print(region)
    # Get the image data array
    data_array = data.get_fdata()
    
    # Initialize combined_data with the same shape as data_array
    if combined_data is None:
        combined_data = np.zeros_like(data_array)
    
    # Create a mask for the current region
    region_mask = data_array > 0
    
    # Replace the values in the region with the corresponding expansion factor
    combined_data[region_mask] = expansion_factors[count]
    
    # Create a NIfTI image from the modified data array for visualization
    img = nib.Nifti1Image(region_mask.astype(np.float32) * expansion_factors[count], affine=data.affine, header=data.header)
    img_bin = binarize_img(img)
    # plt.figure(figsize=(10, 5))
    # plot_roi(img_bin, draw_cross=False, cmap=cmap[count], alpha=1)
    # plt.savefig('./rois/regionROI_' + region + '.png', dpi=1200)
    # plt.savefig('./rois/regionROI_' + region + '.svg', dpi=1200)
    # plt.show()
    count += 1

# Create a new NIfTI image from the combined data
combined_img = nib.Nifti1Image(combined_data, affine=data.affine, header=data.header)

# Save the combined NIfTI image
nib.save(combined_img, './rois/combined_expansion_factors.nii')

print(np.unique(combined_data))

# Plot the combined NIfTI image
plt.figure(figsize=(10, 5))
plot_roi(combined_img, draw_cross=False)
plt.savefig('./rois/combined_expansion_factors.png', dpi=1200)
plt.savefig('./rois/combined_expansion_factors.svg', dpi=1200)
plt.show()
