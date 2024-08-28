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
        'V1': 'occipital', 'V2': 'occipital', 'V3': 'occipital',
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

# Fetch Neurosynth data
neurosynth_db = get_neurosynth_db()

# Example usage
dset = convert_neurosynth_to_dataset(
    coordinates_file=neurosynth_db["coordinates"],
    metadata_file=neurosynth_db["metadata"],
    annotations_files=neurosynth_db["features"],
)
print(dset)

# for region, data in combined_rois.items():
#     output_file = os.path.join('C:/Users/eemeyer/Documents/GitHub/RetinoAreaComparisons/neurosynth/rois', f'{region}_region.png')
#     plot_roi(data, draw_cross=False, output_file=output_file)
#     print(f"Saved combined ROI for {region} as {output_file}")

#     # Save as NIfTI
#     output_nii_file = os.path.join('C:/Users/eemeyer/Documents/GitHub/RetinoAreaComparisons/neurosynth/rois', f'{region}_combined.nii.gz')
#     nib.save(data, output_nii_file)
#     print(f"Saved combined ROI for {region} as {output_nii_file}")

# Create a figure with subplots
regions = list(combined_rois.keys())
fig, axes = plt.subplots(len(regions), 1, figsize=(10, len(regions) * 5))
count = 0

for region, data in combined_rois.items():
    print(region)
    # Get the image data array
    data_array = data.get_fdata()
    
    # # Check if the affine matrices match
    # if not np.allclose(data.affine, dset.masker.mask_img.affine):
    #     print(f"Resampling {region} to match dataset affine.")
    #     data = resample_img(data, target_affine=dset.masker.mask_img.affine)
    #     data_array = data.get_fdata()
    
    # Create a NIfTI image from the chunked data
    img = nib.Nifti1Image(data_array, affine=data.affine, header=data.header)
    img_bin = binarize_img(img)
    # plt.figure(figsize=(10, 5))
    # plot_roi(img, draw_cross=False)
    # plt.show()

    ids = []
    nbin = 100
    for i in range(nbin):
        length_dat = np.round(14371/nbin)
        if i == nbin-1:
            idx = [int(length_dat*i), 14371]
        else:
            idx = [int(length_dat*i), int(length_dat*(i+1))]
        all_ids = dset.ids
        label_slice = all_ids[idx[0]:idx[1]]
        dset_slice = dset.slice(label_slice)

        # Get studies by mask using the NIfTI image
        ids = ids + dset_slice.get_studies_by_mask(img_bin)

    # Run the decoder
    decoder = discrete.NeurosynthDecoder(correction=None)
    decoder.fit(dset)
    decoded_df = decoder.transform(ids=ids)
    print(decoded_df.sort_values(by="probReverse", ascending=False).head())
    tmp = decoded_df.sort_values(by="probReverse", ascending=False)[:20]
    
    index_list = [[idx[22:]] for idx in tmp.index]
    print(index_list)

    # This method decodes the ROI image directly, rather than comparing subsets of the Dataset like the
    # other two.

    # # Load the NIfTI file
    # nii_file_path = os.path.join(nii_dir, region + '_combined.nii.gz')
    # nii_img = nib.load(nii_file_path)

    # # Create a new NIfTI image with the binary mask
    # binary_nii_img = nib.Nifti1Image(nii_img, nii_img.affine, nii_img.header)

    # # Second decoding process
    # logging.info(f"Decoding ROI...")
    # decoder = discrete.ROIAssociationDecoder(img_bin)
    # decoder.fit(dset)

    # # The `transform` method doesn't take any parameters.
    # decoded_df = decoder.transform()

    # print(decoded_df.sort_values(by="r", ascending=False).head())
    # tmp = decoded_df.sort_values(by="r", ascending=False)[:20]
    
    # index_list = [[idx[22:]] for idx in tmp.index]
    # print(index_list)

    # Plot the DataFrame
    axes[count].axis('off')
    axes[count].table(cellText=index_list, colLabels=['Index'], loc='center')
    axes[count].set_title(region)
    count += 1

# Adjust layout and save the figure
plt.tight_layout()
plt.savefig('./test_fig/NeurosynthDecoder_topwords.png')
plt.show()