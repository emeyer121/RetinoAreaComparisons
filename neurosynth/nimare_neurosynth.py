import os
from pprint import pprint
import matplotlib.pyplot as plt
from scipy import stats

import nibabel as nib
import numpy as np
from nilearn.plotting import plot_roi
from nilearn.image import concat_imgs, math_img, binarize_img
from nilearn.image import resample_img
from wordcloud import WordCloud
import matplotlib.colors as mcolors
from collections import defaultdict

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

def is_plural(word1, word2):
    """
    Check if word1 is the plural version of word2.
    """
    if word1.endswith('s') and word1[:-1] == word2:
        return True
    if word2.endswith('s') and word2[:-1] == word1:
        return True
    return False

def extract_top_words(df):

    exclusion_words = [
        'occipital','temporal','parietal','gyrus','sulcus', 'gyri','sulci','lateral','ventral','dorsal',
        'medial','inferior','superior','anterior','posterior','cortex','brain','neuron','v1','primary visual',
        'lingual gyrus','lingual','early visual','visual cortices','visual cortex','cuneus','extrastriate',
        'fusiform gyrus','fusiform gyri','ventral visual','occipital lobe','occipital cortex','parahippocampal cortex',
        'fusiform face','face ffa','fusiform','occipito','occipito temporal',
        'occipitotemporal cortex','occipital temporal','parieto occipital','parahippocampal','occipitotemporal',
        'middle occipital','occipital gyrus','temporal occipital','occipital parietal','anterior hippocampus',
        'lateral occipital','v5','mt','eye field','eye fields','parietal lobules','frontal eye','cortex ppc','ppc',
        'superior parietal','posterior parietal','frontoparietal network','intraparietal sulcus','sulcus ips',
        'intraparietal','parietal occipital','frontoparietal','parietal frontal','parietal network','superior inferior',
        'ips','spl','ffa','parietal cortex','fronto parietal','parietal lobes','inferior superior','inferior occipital',
        'lobule','lobules','parietal cortices','prefrontal parietal','wm','memory wm','wm task','visuo spatial'
    ]

    tmp = df.sort_values(by="probReverse", ascending=False)

    index_list = [[idx[22:]] for idx in tmp.index]
    vals = tmp['probReverse'].values

    # Filter out exclusion words
    filtered_index_list = []
    filtered_vals = []
    for idx, val in zip(index_list, vals):
        if idx[0] not in exclusion_words:
            filtered_index_list.append(idx)
            filtered_vals.append(val)

    filtered_words_vals = {}
    
    for word, val in zip(filtered_index_list, filtered_vals):
        if isinstance(word, list):
            word = word[0]  # Ensure word is a string
        if word not in filtered_words_vals:
            filtered_words_vals[word] = val
        else:
            if val > filtered_words_vals[word]:
                filtered_words_vals[word] = val


    # # Check for plural versions and exclude the one with the lower value
    # final_words_vals = {}
    # words_list = list(filtered_words_vals.keys())
    
    # for i, word1 in enumerate(words_list):
    #     if word1 in final_words_vals:
    #         continue
    #     for j, word2 in enumerate(words_list):
    #         if i != j and is_plural(word1, word2):
    #             if filtered_words_vals[word1] >= filtered_words_vals[word2]:
    #                 final_words_vals[word1] = filtered_words_vals[word1]
    #             else:
    #                 final_words_vals[word2] = filtered_words_vals[word2]
    #             break
    #     else:
    #         final_words_vals[word1] = filtered_words_vals[word1]

    filtered_index_list = list(filtered_words_vals.keys())
    filtered_vals = list(filtered_words_vals.values())

    return filtered_index_list, filtered_vals

def process_region(region, data, dset):
    """
    Process a single region and plot the results.
    
    Parameters:
    region (str): Name of the region.
    data (Nifti1Image): NIfTI image data for the region.
    dset (Dataset): Neurosynth dataset.
    axes (Axes): Matplotlib axes for plotting.
    count (int): Index of the subplot.
    
    Returns:
    int: Updated count.
    """
    print(region)
    data_array = data.get_fdata()
    
    # if not np.allclose(data.affine, dset.masker.mask_img.affine):
    #     print(f"Resampling {region} to match dataset affine.")
    #     data = resample_img(data, target_affine=dset.masker.mask_img.affine)
    #     data_array = data.get_fdata()
    
    img = nib.Nifti1Image(data_array, affine=data.affine, header=data.header)
    img_bin = binarize_img(img)

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

    decoder = discrete.NeurosynthDecoder(correction=None)
    decoder.fit(dset)
    decoded_df = decoder.transform(ids=ids)
    print(decoded_df.sort_values(by="probReverse", ascending=False)[:20])
    
    filtered_index_list, filtered_vals = extract_top_words(decoded_df)

    return filtered_index_list, filtered_vals

def generate_wordcloud(words,vals,cmap):
    """
    Generate word clouds for each coordinate set based on TF-IDF scores.
    Displays the word clouds.
    """

    max_font_size = 100
    colormap = plt.get_cmap(cmap)

    def color_func(word, font_size, position, orientation, random_state=None, **kwargs):
        base_color = colormap(0.5)  # Get the base color from the colormap
        lightness = int(font_size / max_font_size * 100)
        r, g, b, _ = base_color
        h, l, s = mcolors.rgb_to_hsv((r, g, b))
        l = lightness / 100.0
        r, g, b = mcolors.hsv_to_rgb((h, l, s))
        return f"rgb({int(r*255)}, {int(g*255)}, {int(b*255)})"
    
    # Create a dictionary from words and vals
    wordcloud_data = {word: val for word, val in zip(words, vals)}

    wordcloud = WordCloud(width=600, height=800, background_color='white', relative_scaling=1,
    prefer_horizontal=1, max_words=20, max_font_size=max_font_size, color_func=color_func).generate_from_frequencies(wordcloud_data)

    return wordcloud

def merge_and_expand_words_vals(all_words, all_vals, expansion_factors):
    """
    Merge words and values from multiple regions, applying expansion factors and taking the top value for each word.
    
    Parameters:
    all_words (dict): Dictionary with region names as keys and lists of words as values.
    all_vals (dict): Dictionary with region names as keys and lists of values as values.
    expansion_factors (dict): Dictionary with region names as keys and expansion factors as values.
    
    Returns:
    dict: Dictionary with words as keys and top values as values.
    """
    merged_dict = defaultdict(float)

    for idx, region in enumerate(all_words.keys()):
        words = all_words[region]
        vals = all_vals[region]
        # print(np.min(vals), np.max(vals))
        # vals = vals / np.max(vals)  # Normalize values
        factor = expansion_factors[idx]
        
        # Consider only the top 50 words and their values
        top_words_vals = sorted(zip(words, vals), key=lambda x: x[1], reverse=True)[:50]
        # Extract just the values from the sorted dictionary
        # values = [value for _, value in top_words_vals]
        # print(values)
        # new_vals = stats.zscore(values)
        # print(new_vals)
        
        for word, val in top_words_vals:
            print(word, val)
            expanded_val = val * factor
            if expanded_val > merged_dict[word]:
                merged_dict[word] = expanded_val

    return dict(merged_dict)
    

# Extract region ROIS
nii_dir = os.path.join(os.path.dirname(__file__), 'rois')
combined_rois = combine_rois(nii_dir)

# Fetch Neurosynth data
neurosynth_db = get_neurosynth_db()

# Convert Neurosynth data to dataset
dset = convert_neurosynth_to_dataset(
    coordinates_file=neurosynth_db["coordinates"],
    metadata_file=neurosynth_db["metadata"],
    annotations_files=neurosynth_db["features"],
)

# Create a figure with subplots
regions = list(combined_rois.keys())
fig, axes = plt.subplots(1, len(regions), figsize=(len(regions) * 5, 10))
all_words = {}
all_vals = {}
count = 0

for region, data in combined_rois.items():
    filtered_list, filtered_vals = process_region(region, data, dset)
    print(filtered_list[:20])
    all_words[region] = filtered_list
    all_vals[region] = filtered_vals

    # Prepare data for the table
    table_data = [[word] for word in filtered_list[:20]]
    col_labels = [region]

    # Plot the DataFrame
    axes[count].axis('off')
    axes[count].table(cellText=table_data, colLabels=col_labels, loc='center')
    count += 1

# Adjust layout and save the figure
plt.tight_layout()
# plt.savefig('./test_fig/BrainmapDecoder_topwords_filteredv4.png',dpi=1200)
# plt.savefig('./test_fig/BrainmapDecoder_topwords_filteredv4.svg',dpi=1200)
plt.show()


cmap = ['blue','firebrick','orange','green','purple']

# Create a figure with subplots for histograms
fig_hist, axes_hist = plt.subplots(1, len(regions), figsize=(len(regions) * 5, 10))
count = 0

for idx, region in enumerate(regions):
    vals = all_vals[region]
    cutoff = vals[19] if len(vals) > 19 else vals[-1]  # Ensure there are at least 20 values

    # Plot histogram
    axes_hist[count].hist(vals, bins=30, color='gray', alpha=0.7)
    axes_hist[count].axvline(cutoff, color=cmap[idx], linestyle='dashed', linewidth=2)
    axes_hist[count].set_title(region)
    axes_hist[count].set_xlabel('Values')
    axes_hist[count].set_ylabel('Frequency')
    count += 1

# Adjust layout and display the histogram plot
plt.tight_layout()
plt.show()

cmap = ['Blues','Reds','Oranges','Greens','BuGn_r']
count = 0

fig, axes = plt.subplots(1, len(regions), figsize=(20, 10))
for region, data in combined_rois.items():
    wordcloud = generate_wordcloud(all_words[region],all_vals[region],cmap[count])

    axes[count].imshow(wordcloud, interpolation='bilinear')
    axes[count].set_title(f'Word Cloud for {region}')
    axes[count].axis('off')

    count += 1

# Adjust layout and save the figure
plt.tight_layout()
# plt.savefig('./test_fig/BrainmapDecoder_regioncloudv5.png',dpi=1200)
# plt.savefig('./test_fig/BrainmapDecoder_regioncloudv5.svg',dpi=1200)
plt.show()

expansion_factors = [1.3, 2.67, 3.94, 4.25, 7.11]

merged_words_vals = merge_and_expand_words_vals(all_words, all_vals, expansion_factors)

# Sort the dictionary by values in descending order and get the top 50 items
sorted_merged_words_vals = sorted(merged_words_vals.items(), key=lambda item: item[1], reverse=True)
table_data = sorted_merged_words_vals[:50]
# print(table_data)

# Create a new figure and axis for the table
fig, ax = plt.subplots()
ax.axis('off')  # Hide the axis

# Create the table
table = ax.table(cellText=table_data, colLabels=['Word', 'Value'], loc='center')

# Adjust layout and save the figure
plt.tight_layout()
# plt.savefig('./test_fig/BrainmapDecoder_topwords_expansion_norm_v4.png', dpi=1200)
plt.show()

# Initialize a dictionary to count the number of regions each word appears in the top 20
term_region_counts = {word: 0 for word in merged_words_vals.keys()}

# Iterate over each region to find the top 20 words and update the counts
for region in all_words.keys():
    top_20_words_in_region = [word for word, val in sorted(zip(all_words[region], all_vals[region]), key=lambda x: x[1], reverse=True)[:20]]
    for word in merged_words_vals.keys():
        if word in top_20_words_in_region:
            term_region_counts[word] += 1

# Define a function to map the count to a color
def term_color_func(word, *args, **kwargs):
    count = term_region_counts.get(word, 0)
    # Map count to a shade of gray (0: black, 5: light gray)
    gray_value = int(255 * (count / 5))
    return f'rgb({gray_value}, {gray_value}, {gray_value})'

wordcloud_exp = WordCloud(width=1200, height=600, background_color='white',
    prefer_horizontal=1, max_words=50, color_func=term_color_func).generate_from_frequencies(merged_words_vals)

plt.figure(figsize=(10, 5))
plt.imshow(wordcloud_exp, interpolation='bilinear')
plt.axis('off')
plt.title('Terms Associated with Cortical Expansion')
# plt.tight_layout(pad=0)
# plt.savefig('./test_fig/BrainmapDecoder_expansioncloud_norm_v4.png',dpi=1200)
# plt.savefig('./test_fig/BrainmapDecoder_expansioncloud_norm_v4.svg',dpi=1200)
plt.show()


# new_word_dict = {'visual':0.351, 'spatial':0.267,'vision':0.231,'saccades':0.218,'attentional':0.213,
#                  'eye movements':0.188,'attention':0.187,'visuospatial':0.187,'rotation':0.179,'sighted':0.172,
#                  'spatial attention':0.167,'visuo':0.163,'motion':0.157,'location':0.146,'calculation':0.146,
#                  'saccade':0.142,'objects':0.139,'object':0.137,'selective attention':0.127,'fixation':0.121,
#                  'target':0.108,'color':0.106,'spatial information':0.095,'task':0.094,'visual motion':0.089,
#                  'tasks':0.088,'visual attention':0.088,'arithmetic':0.083,'locations':0.083,'dimensional':0.082,
#                  'orienting':0.08,'attention network':0.076,'attended':0.075,'shifts':0.072,'navigation':0.069,
#                  'depth':0.068,'video':0.067,'attentional control':0.059,'videos':0.059,'load':0.058,'visually presented':0.058,
#                  'imagery':0.057,'visual stream':0.057,'grasping':0.056,'preparatory':0.055,'action observation':0.054,
#                  'virtual':0.053,'memory load':0.052,'space':0.051,'subtraction':0.051,'working memory':0.051,'working':0.05,
#                  'body':0.049,'gaze':0.048,'mental imagery':0.048,'performance':0.047,'passively':0.045,'rule':0.044,
#                  'orientation':0.044,'perceptual':0.043,'action':0.042,'multisensory':0.042,'actions':0.041,'spatially':0.04}

# wordcloud_exp = WordCloud(width=1000, height=500, background_color='white',
#     prefer_horizontal=1, max_words=50, relative_scaling=1, color_func=term_color_func).generate_from_frequencies(new_word_dict)

# plt.figure(figsize=(10, 5))
# plt.imshow(wordcloud_exp, interpolation='bilinear')
# plt.axis('off')
# plt.title('Terms Associated with Cortical Expansion')
# plt.tight_layout(pad=0)
# plt.savefig('./test_fig/NeuroVaultDecoder_expansioncloud_v2.png',dpi=1200)
# plt.savefig('./test_fig/NeuroVaultDecoder_expansioncloud_v2.svg',dpi=1200)
# plt.show()
