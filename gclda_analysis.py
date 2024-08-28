import numpy as np
import pandas as pd
import random
import logging
from wordcloud import WordCloud
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import inflect


iter = '1000'

# Load activation assignments
#activation_assignments = pd.read_csv('C:/Users/eemeyer/Documents/GitHub/RetinoAreaComparisons/python_gclda/examples/gclda_results/2015Filtered2_TrnTst1p1_100T_2R_alpha0.100_beta0.010_gamma0.010_delta1.000_25dobs_50.0roi_1symmetric_1/iter_'+iter+'/ActivationAssignments.csv')

# Load topicXword count matrix
#topic_word_count = pd.read_csv('C:/Users/eemeyer/Documents/GitHub/RetinoAreaComparisons/python_gclda/examples/gclda_res#ults/2015Filtered2_TrnTst1p1_100T_2R_alpha0.100_beta0.010_gamma0.010_delta1.000_25dobs_50.0roi_1symmetric_1/iter_'+iter+'/Topic_X_Word_CountMatrix.csv')


# Load activation assignments
activation_assignments = pd.read_csv('C:/Users/eemeyer/Documents/GitHub/RetinoAreaComparisons/python_gclda/example_results/2015Filtered2_TrnTst1P1_100T_2R_alpha0.100_beta0.010_gamma0.010_delta1.000_25dobs_50.0roi_1symmetric_1/iter_'+iter+'/ActivationAssignments.csv')

# Load topicXword count matrix
topic_word_count = pd.read_csv('C:/Users/eemeyer/Documents/GitHub/RetinoAreaComparisons/python_gclda/example_results/2015Filtered2_TrnTst1P1_100T_2R_alpha0.100_beta0.010_gamma0.010_delta1.000_25dobs_50.0roi_1symmetric_1/iter_'+iter+'/Topic_X_Word_CountMatrix.csv')

print(activation_assignments.head())
print(topic_word_count.head())

def define_coordinate_sets(file_path="neurosynth\VTPMcoords.txt"):
    """
    Define coordinate sets for different brain regions and their scaling factors.
    Returns a dictionary with region names, coordinates, and scaling factors.
    """
    occipitalROI = list(range(2, 8)) + [14]
    dorsalROI = list(range(8, 10))
    lateralROI = list(range(10, 14))
    ventralROI = list(range(15, 19))
    parietalROI = list(range(19, 26))

    roi_groups = {
        'occipitalROI': [],
        'dorsalROI': [],
        'lateralROI': [],
        'ventralROI': [],
        'parietalROI': []
    }
    coordinate_sets = {}
    try:
        with open(file_path, 'r') as file:
            for count, line in enumerate(file):
                parts = line.strip().split()
                coordinates = [(float(parts[0]), float(parts[1]), float(parts[2]))]
                ROI = float(parts[3])
                coordinate_sets[count] = {"coordinates": coordinates, "ROI": ROI}
    except Exception as e:
        logging.error(f"Error reading coordinate sets from file: {e}")
        raise

    for data in coordinate_sets.items():
        if data[1]['ROI'] in occipitalROI:
            roi_groups['occipitalROI'].append(data[1]['coordinates'][0])
        elif data[1]['ROI'] in dorsalROI:
            roi_groups['dorsalROI'].append(data[1]['coordinates'][0])
        elif data[1]['ROI'] in lateralROI:
            roi_groups['lateralROI'].append(data[1]['coordinates'][0])
        elif data[1]['ROI'] in ventralROI:
            roi_groups['ventralROI'].append(data[1]['coordinates'][0])
        elif data[1]['ROI'] in parietalROI:
            roi_groups['parietalROI'].append(data[1]['coordinates'][0])
    
    for key in roi_groups:
        print(len(roi_groups[key]))
        print(roi_groups[key][1:10])
    #print(roi_groups['ventralROI'])

    return {
        "occipital": {"coordinates": roi_groups['occipitalROI'], "scaling_factor": 1.3},
        "ventral": {"coordinates": roi_groups['ventralROI'], "scaling_factor": 2.67},
        "dorsal": {"coordinates": roi_groups['dorsalROI'], "scaling_factor": 3.94},
        "lateral": {"coordinates": roi_groups['lateralROI'], "scaling_factor": 4.25},
        "parietal": {"coordinates": roi_groups['parietalROI'], "scaling_factor": 7.11}
    }


def match_coordinates(activation_assignments, coordinate_sets):
    """
    Matches the coordinates to the x, y, and z coordinates in activation assignments
    and stores the topics associated with each region.
    Returns a dictionary with region names and their associated topics.
    """
    region_topics = {}
    for region, data in coordinate_sets.items():
        coordinates = data["coordinates"]
        topics = []
        for coord in coordinates:
            x, y, z = coord
            # Find the matching coordinates in activation assignments
            matching_coords = activation_assignments[
                (activation_assignments.iloc[:, 0] == x) &
                (activation_assignments.iloc[:, 1] == y) &
                (activation_assignments.iloc[:, 2] == z)
            ]
            if not matching_coords.empty: #and np.shape(matching_coords)[0] == 1:
                # Get the topics associated with the matching coordinates
                matching_topics = matching_coords.iloc[:, 3].values
                topics.extend(matching_topics)
        region_topics[region] = topics
    return region_topics

def calculate_word_frequencies(topic_word_count, region_topics):
    """
    Calculates the frequency of words associated with each topic in each region.
    Returns a dictionary with region names and their associated word frequencies.
    """

    region_word_frequencies = {}
    for region, topics in region_topics.items():
        word_frequencies = np.zeros(topic_word_count.shape[0])
        for topic in topics:
            column_name = f"Topic_{str(topic).zfill(2)}"
            if column_name in topic_word_count.columns:
                word_frequencies += topic_word_count[column_name].values
            else:
                print(f"Warning: Column {column_name} not found in topic_word_count")
        region_word_frequencies[region] = word_frequencies
    return region_word_frequencies

def remove_common_top_words(region_word_frequencies, topic_word_count):
    """
    Removes words that have the top 20 frequency in all regions iteratively until there are no words that appear in all five regions.
    """

    words = topic_word_count.iloc[:, 0].values
    p = inflect.engine()
    
    while True:
        # Find the top 20 words in each region
        top_words_in_regions = []
        for region, frequencies in region_word_frequencies.items():
            frequencies_series = pd.Series(frequencies, index=words)
            top_words_in_region = set(frequencies_series.nlargest(20).index)
            top_words_in_regions.append(top_words_in_region)

        # Find common top words in all regions
        common_top_words = set.intersection(*top_words_in_regions)
        if not common_top_words:
            break

        # Remove common top words from frequencies
        for region in region_word_frequencies:
            for word in common_top_words:
                word_index = np.where(words == word)[0][0]
                region_word_frequencies[region][word_index] = 0

    return region_word_frequencies

def create_word_clouds(word_frequencies, topic_word_count):
    """
    Creates a word cloud for each region using the word frequencies and the words from the first column of topic_word_count.
    """

    cmap = ['Blues','Reds','Oranges','Greens','BuGn_r']
    logging.info("Generating word clouds...")
    fig, axes = plt.subplots(1, len(coordinate_sets), figsize=(20, 10))

    words = topic_word_count.iloc[:, 0].values  # Get the words from the first column

    idx = 0
    for region, frequencies in word_frequencies.items():
        # Create a dictionary of words and their corresponding frequencies
        word_freq_dict = {word: freq for word, freq in zip(words, frequencies)}

        max_font_size = 300
        colormap = plt.get_cmap(cmap[idx])

        def color_func(word, font_size, position, orientation, random_state=None, **kwargs):
            base_color = colormap(0.5)  # Get the base color from the colormap
            lightness = int(font_size / max_font_size * 100)
            r, g, b, _ = base_color
            h, l, s = mcolors.rgb_to_hsv((r, g, b))
            l = lightness / 100.0
            r, g, b = mcolors.hsv_to_rgb((h, l, s))
            return f"rgb({int(r*255)}, {int(g*255)}, {int(b*255)})"

        wordcloud = WordCloud(width=600, height=800, background_color='white',
        prefer_horizontal=1, max_words=25, color_func=color_func).generate_from_frequencies(word_freq_dict)
        axes[idx].imshow(wordcloud, interpolation='bilinear')
        axes[idx].set_title(f'Word Cloud for {region}')
        axes[idx].axis('off')
        idx += 1
    plt.tight_layout()
    plt.show()

def calculate_expansion_weighted_scores(word_frequencies, topic_word_count, coordinate_sets):
    """
    Calculate expansion-weighted scores for each term.
    Returns a Series of weighted scores.
    """

    words = topic_word_count.iloc[:, 0].values  # Get the words from the first column

    # Initialize a dictionary to accumulate word frequencies across all regions
    accumulated_word_freq_dict = {word: 0 for word in words}

    for region, frequencies in word_frequencies.items():
        # Accumulate the frequencies for each word
        for word, freq in zip(words, frequencies):
            accumulated_word_freq_dict[word] += freq

    expansion_factors = {region: data['scaling_factor'] for region, data in coordinate_sets.items()}
    weighted_scores = pd.Series(accumulated_word_freq_dict)

    for region, frequencies in word_frequencies.items():
        # Normalize the frequencies for the region
        total_freq = sum(frequencies)
        normalized_frequencies = [freq / total_freq for freq in frequencies]

        for word, norm_freq in zip(words, normalized_frequencies):
            weighted_scores[word] += norm_freq * expansion_factors[region]

    return weighted_scores

def generate_expansion_word_cloud(word_frequencies, topic_word_count, coordinate_sets):
    """
    Generate a word cloud of terms associated with cortical expansion.
    Displays the word cloud and logs top terms.
    """
    logging.info("Generating expansion-associated word cloud...")
    
    weighted_scores = calculate_expansion_weighted_scores(word_frequencies, topic_word_count, coordinate_sets)
    
    # Normalize scores to be non-negative
    min_score = weighted_scores.min()
    if min_score < 0:
        weighted_scores = weighted_scores - min_score

     # Initialize a dictionary to count the number of regions each word appears in the top 30
    term_region_counts = {word: 0 for word in weighted_scores.index}
    
    # Iterate over each region to find the top 30 words and update the counts
    for region, frequencies in word_frequencies.items():
        frequencies_series = pd.Series(frequencies, index=topic_word_count.iloc[:, 0].values)
        top_50_words_in_region = frequencies_series.nlargest(50).index
        for word in top_50_words_in_region:
            if word in term_region_counts:
                term_region_counts[word] += 1
        
    # Define a function to map the count to a color
    def term_color_func(word, *args, **kwargs):
        count = term_region_counts.get(word, 0)
        # Map count to a shade of gray (0: black, 5: light gray)
        gray_value = int(255 * ((count-1) / 5))
        return f'rgb({gray_value}, {gray_value}, {gray_value})'
    
    # Create and display the word cloud
    wordcloud = WordCloud(width=800, height=400, background_color='white', max_words=50, color_func=term_color_func).generate_from_frequencies(weighted_scores)
    
    plt.figure(figsize=(10, 5))
    plt.imshow(wordcloud, interpolation='bilinear')
    plt.axis('off')
    plt.title('Terms Associated with Cortical Expansion')
    plt.tight_layout(pad=0)
    plt.show()
    
    # Log top terms
    top_terms = weighted_scores.nlargest(20)
    logging.info("Top 20 terms associated with cortical expansion:")
    for term, score in top_terms.items():
        logging.info(f"{term}: {score:.4f}")


coordinate_sets = define_coordinate_sets()

# Call the function to match coordinates
region_topics = match_coordinates(activation_assignments, coordinate_sets)

# Call the function to calculate word frequencies
word_frequencies = calculate_word_frequencies(topic_word_count, region_topics)

filtered_word_frequencies = remove_common_top_words(word_frequencies, topic_word_count)

# Print the word frequencies for each region
for region, frequencies in filtered_word_frequencies.items():
    print(f"Region: {region}, # of topics: {len(region_topics[region])}")
    #print(region_topics[region])
    print(f"Word Frequencies: {frequencies}")
    print()

# Create word clouds for each region
create_word_clouds(filtered_word_frequencies, topic_word_count)

generate_expansion_word_cloud(filtered_word_frequencies, topic_word_count, coordinate_sets)
