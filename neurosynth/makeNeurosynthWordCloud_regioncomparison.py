#pip install numpy pandas scipy scikit-learn wordcloud matplotlib seaborn plotly neurosynth

# TF-IDF stands for Term Frequency-Inverse Document Frequency. 
# It's a statistical measure used to evaluate the importance of a word in a document within a collection or 
# corpus of documents. In the context of your neuroimaging analysis, each brain region is treated as a "document," 
# and the terms associated with that region are the "words" in that document.
# Here's a breakdown of TF-IDF:
# Term Frequency (TF):
#
# This is simply the number of times a term appears in a document (in this case, a brain region), 
# normalized by the total number of terms in that document.
# TF(t,d) = (Number of times term t appears in document d) / (Total number of terms in document d)
#
#
# Inverse Document Frequency (IDF):
# This measures how common or rare a term is across all documents (brain regions).
# IDF(t) = log(Total number of documents / Number of documents containing term t)
# TF-IDF:
# The TF-IDF score is the product of TF and IDF.
# TF-IDF(t,d) = TF(t,d) * IDF(t)
#
# The key ideas behind TF-IDF are:

# Terms that are common in a particular document (brain region), but rare across all documents, 
# receive a higher TF-IDF score.
# Terms that are common across all documents receive a lower score.

# In your neuroimaging context:

# A high TF-IDF score for a term in a specific brain region suggests that this term is particularly important 
# or characteristic for that region.
# Terms with high TF-IDF scores are likely to be more discriminative and informative about the specific functions 
# or properties of that brain region.
# Terms with low TF-IDF scores are likely to be common across many brain regions and thus less informative about 
# any specific region.
# By using TF-IDF, you're able to identify terms that are uniquely or strongly associated with each brain region, 
# which can provide insights into the functional specialization of different areas of the brain. 
# This approach helps to filter out common, less informative terms and highlight those that are more likely to be 
# meaningful in characterizing each region's specific role or properties.


import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.feature_extraction.text import TfidfTransformer
from wordcloud import WordCloud
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency 
from scipy.stats import wilcoxon
import plotly.graph_objects as go
import logging
import random

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_data():
    """
    Load Neurosynth data from files.
    Returns metadata, coordinates, and feature DataFrames.
    """
    logging.info("Loading data...")
    try:
        feature_data_sparse = sparse.load_npz(".\data-neurosynth_version-7_vocab-terms_source-abstract_type-tfidf_features.npz")
        feature_data = feature_data_sparse.todense()
        metadata_df = pd.read_table(".\data-neurosynth_version-7_metadata.tsv.gz")
        coordinates_df = pd.read_table(".\data-neurosynth_version-7_coordinates.tsv.gz")
        feature_names = np.genfromtxt(".\data-neurosynth_version-7_vocab-terms_vocabulary.txt", dtype=str, delimiter="\t").tolist()
        
        ids = metadata_df["id"].tolist()
        feature_df = pd.DataFrame(index=ids, columns=feature_names, data=feature_data)
        
        return metadata_df, coordinates_df, feature_df
    except Exception as e:
        logging.error(f"Error loading data: {e}")
        raise

def define_coordinate_sets(file_path=".\VTPMcoords.txt"):
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

    new_array = []
    for data in coordinate_sets.items():
        if data[1]['ROI'] in occipitalROI:
            roi_groups['occipitalROI'].append(data[1]['coordinates'][0])
            new_array.append([data[1]['coordinates'][0][0], data[1]['coordinates'][0][1], data[1]['coordinates'][0][2], data[1]['ROI'], 1.3])
        elif data[1]['ROI'] in dorsalROI:
            roi_groups['dorsalROI'].append(data[1]['coordinates'][0])
            new_array.append([data[1]['coordinates'][0][0], data[1]['coordinates'][0][1], data[1]['coordinates'][0][2], data[1]['ROI'], 2.67])
        elif data[1]['ROI'] in lateralROI:
            roi_groups['lateralROI'].append(data[1]['coordinates'][0])
            new_array.append([data[1]['coordinates'][0][0], data[1]['coordinates'][0][1], data[1]['coordinates'][0][2], data[1]['ROI'], 3.94])
        elif data[1]['ROI'] in ventralROI:
            roi_groups['ventralROI'].append(data[1]['coordinates'][0])
            new_array.append([data[1]['coordinates'][0][0], data[1]['coordinates'][0][1], data[1]['coordinates'][0][2], data[1]['ROI'], 4.25])
        elif data[1]['ROI'] in parietalROI:
            roi_groups['parietalROI'].append(data[1]['coordinates'][0])
            new_array.append([data[1]['coordinates'][0][0], data[1]['coordinates'][0][1], data[1]['coordinates'][0][2], data[1]['ROI'], 7.11])

    coordinate_sets = {
        "occipital": {"coordinates": [roi_groups['occipitalROI'][i] for i in random.sample(range(len(roi_groups['occipitalROI'])), 100)], "scaling_factor": 1.3},
        "ventral": {"coordinates": [roi_groups['ventralROI'][i] for i in random.sample(range(len(roi_groups['ventralROI'])), 100)], "scaling_factor": 2.67},
        "dorsal": {"coordinates": [roi_groups['dorsalROI'][i] for i in random.sample(range(len(roi_groups['dorsalROI'])), 100)], "scaling_factor": 3.94},
        "lateral": {"coordinates": [roi_groups['lateralROI'][i] for i in random.sample(range(len(roi_groups['lateralROI'])), 100)], "scaling_factor": 4.25},
        "parietal": {"coordinates": [roi_groups['parietalROI'][i] for i in random.sample(range(len(roi_groups['parietalROI'])), 100)], "scaling_factor": 7.11}
    }

    return coordinate_sets, new_array

#    print(roi_groups)
#    print(roi_groups['occipitalROI'])

#    return coordinate_sets

#    return {
#        "occipital": {"coordinates": [(-10, -94, -2), (14, -94, -2)], "scaling_factor": 1.5},
#        "lateral": {"coordinates": [(44, -66, 0), (-46, -72, 0)], "scaling_factor": 4},
#        "parietal": {"coordinates": [(-32, -54, 52), (36, -46, 52)], "scaling_factor": 8},
#        "ventral": {"coordinates": [(-38, -50, -20), (48, -50, -24)], "scaling_factor": 2.5}
#    }

def extract_terms_z(coordinates_df, feature_df, coordinates, scaling_factor, radius=6):
    """
    Extract terms for a given set of coordinates within a specified radius.
    Returns a Series of term scores.
    """
    closest_studies = []
    for idx, coord in enumerate(coordinates):
        print(idx, coord)
        x, y, z = coord
        distances = coordinates_df.apply(lambda row: np.sqrt((row['x'] - x)**2 + (row['y'] - y)**2 + (row['z'] - z)**2), axis=1)
        studies_within_radius = coordinates_df[distances <= radius]
        closest_studies.extend(studies_within_radius['id'].tolist())
    terms_z = feature_df.loc[closest_studies].mean(axis=0) #* scaling_factor
    logging.info(f"Ending extraction")
    return terms_z

def merge_singular_plural(terms):
    """
    Merge singular and plural forms of terms.
    Returns a Series with merged term scores.
    """
    singular_plural_map = {}
    for term in terms.index:
        if term.endswith('s'):
            singular_term = term[:-1]
            if singular_term in terms.index:
                singular_plural_map[term] = singular_term
            else:
                singular_plural_map[term] = term
        else:
            plural_term = term + 's'
            if plural_term in terms.index:
                singular_plural_map[plural_term] = term
            singular_plural_map[term] = term
    # HAD TO ADD IN DTYPE OTHERWISE IT MADE VALUES INT AND ROUNDED EVERYTHING TO 0
    merged_terms = pd.Series(0, index=set(singular_plural_map.values()), dtype=float)
    for term, merged_term in singular_plural_map.items():
        merged_terms[merged_term] += terms[term]
    return merged_terms

def process_terms(metadata_df, coordinates_df, feature_df, coordinate_sets):
    """
    Process terms for each coordinate set, excluding specific words and merging singular/plural forms.
    Returns a DataFrame of processed term data for all coordinate sets.
    """
    logging.info("Processing terms...")
    exclusion_words = [
        'occipital', 'temporal', 'parietal', 'gyrus', 'sulcus', 'gyri', 'sulci', 'lateral',
        'ventral', 'dorsal', 'medial', 'inferior', 'superior', 'anterior', 'posterior', 
        'cortex', 'brain', 'visual', 'network', 'connectivity', 'task', 'tasks', 'ppc', 'frontal',
        'functional magnetic','signal','adult','resonance','stimulus','stimuli','magnetic resonance','using',
        'magnetic','using','involved','state','human','individual','information','level','cue',
        'performance','condition','age','pattern','involved','suggest','response','control','effect','group',
        'groups','adults','function','time','role','greater','model','trial','trials','activation','regions',
        'mechanism','mechanisms','study','studies','results','result','subjects','subject','processing',
        'conditions','cortex','cortical','neural','neurons','neuron','neuroimaging','neuroscience','effects',
        'processes','representation','activations','self','non'
    ]

    term_data_z = {}
    top_terms = {}
    for set_name, data in coordinate_sets.items():
        terms_z = extract_terms_z(coordinates_df, feature_df, data["coordinates"], data["scaling_factor"])
        terms_z = terms_z.drop(exclusion_words, errors='ignore')
        terms_z = merge_singular_plural(terms_z)
        term_data_z[set_name] = terms_z
        top_terms[set_name] = terms_z.nlargest(20).index.tolist()
        logging.info(f"Extracted {len(terms_z)} terms for {set_name}")
        logging.info(f"Top 10 terms for {set_name}: {terms_z.nlargest(10).index.tolist()}")

    # Count occurrences of each term across all sets
    term_counts = {}
    for terms in top_terms.values():
        for term in terms:
            if term not in term_counts:
                term_counts[term] = 0
            term_counts[term] += 1

    # Identify terms that appear in at least 4 different sets
    common_terms = {term for term, count in term_counts.items() if count >= 4}
    print(common_terms)

    # Remove common terms from each set
    for set_name in term_data_z:
        term_data_z[set_name] = term_data_z[set_name].drop(common_terms, errors='ignore')
        logging.info(f"Top 10 terms for {set_name}: {term_data_z[set_name].nlargest(10).index.tolist()}")

    for set_name in term_data_z:
        top_terms[set_name] = term_data_z[set_name].nlargest(20).index.tolist()

    term_counts = {}
    for terms in top_terms.values():
        for term in terms:
            if term not in term_counts:
                term_counts[term] = 0
            term_counts[term] += 1
    
    # Identify terms that appear in at least 4 different sets
    common_terms = {term for term, count in term_counts.items() if count >= 4}
    print(common_terms)

    # Remove common terms from each set
    for set_name in term_data_z:
        term_data_z[set_name] = term_data_z[set_name].drop(common_terms, errors='ignore')
        logging.info(f"Top 10 terms for {set_name}: {term_data_z[set_name].nlargest(10).index.tolist()}")


    combined_term_data = pd.DataFrame(term_data_z).fillna(0)
    logging.info(f"Combined term data shape: {combined_term_data.shape}")
    return combined_term_data

def calculate_tfidf(combined_term_data):
    """
    Calculate TF-IDF scores for the combined term data.
    Returns a DataFrame of TF-IDF scores.
    """
    logging.info("Calculating TF-IDF scores...")
    tfidf_transformer = TfidfTransformer(norm=None)
    tfidf_scores = tfidf_transformer.fit_transform(combined_term_data)
    return pd.DataFrame(tfidf_scores.toarray(), index=combined_term_data.index, columns=combined_term_data.columns)

def generate_word_clouds(tfidf_df, coordinate_sets):
    """
    Generate word clouds for each coordinate set based on TF-IDF scores.
    Displays the word clouds.
    """
    cmap = ['Blues','Reds','Oranges','Greens','BuGn_r']
    logging.info("Generating word clouds...")
    fig, axes = plt.subplots(1, len(coordinate_sets), figsize=(20, 10))

    for idx, set_name in enumerate(coordinate_sets.keys()):
        if set_name not in tfidf_df.columns:
            logging.warning(f"No data for {set_name}, skipping word cloud generation")
            axes[idx].text(0.5, 0.5, f"No data for {set_name}", ha='center', va='center')
            axes[idx].axis('off')
            continue
        tfidf_scores = tfidf_df[set_name]
        wordcloud_data = {term: score for term, score in tfidf_scores.items() if score > 0}
        if not wordcloud_data:
            logging.warning(f"No positive scores for {set_name}, skipping word cloud generation")
            axes[idx].text(0.5, 0.5, f"No positive scores for {set_name}", ha='center', va='center')
            axes[idx].axis('off')
            continue

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
        prefer_horizontal=1, max_words=20, color_func=color_func).generate_from_frequencies(wordcloud_data)
        axes[idx].imshow(wordcloud, interpolation='bilinear')
        axes[idx].set_title(f'Word Cloud for {set_name}')
        axes[idx].axis('off')
    plt.tight_layout()
    plt.show()

def generate_heatmap(tfidf_df):
    """
    Generate a heatmap of TF-IDF scores for top terms across coordinate sets.
    Displays the heatmap.
    """
    logging.info("Generating heatmap...")
    top_terms = tfidf_df.mean(axis=1).sort_values(ascending=False).head(20).index
    filtered_tfidf_df = tfidf_df.loc[top_terms]
    
    plt.figure(figsize=(14, 10))
    sns.heatmap(filtered_tfidf_df, annot=True, cmap="coolwarm", cbar=True, linewidths=0.5)
    plt.title('Heatmap of TF-IDF Scores Across Coordinate Sets for Top Terms')
    plt.xlabel('Coordinate Sets')
    plt.ylabel('Terms')
    plt.show()

def generate_bar_graph(tfidf_df):
    """
    Generate a bar graph of top terms by mean TF-IDF score across coordinate sets.
    Displays the bar graph and logs the top terms.
    """
    logging.info("Generating bar graph...")
    
    # Calculate mean TF-IDF score across coordinates for each term
    mean_tfidf_scores = tfidf_df.mean(axis=1)
    
    # Get top 20 terms by mean TF-IDF score
    top_terms = mean_tfidf_scores.nlargest(20)
    
    plt.figure(figsize=(14, 8))
    top_terms.plot(kind='bar')
    plt.title('Top 20 Terms by Mean TF-IDF Score Across Coordinate Sets')
    plt.xlabel('Terms')
    plt.ylabel('Mean TF-IDF Score')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()
    
    # Log the top terms
    logging.info("Top 20 terms by mean TF-IDF score:")
    for term, score in top_terms.items():
        logging.info(f"  {term}: {score:.4f}")

def analyze_tfidf_distribution(tfidf_df):
    """
    Analyze the distribution of TF-IDF scores across regions.
    Displays a boxplot and logs basic statistics.
    """
    logging.info("Analyzing TF-IDF score distribution...")
    
    # Calculate basic statistics for each region
    stats = tfidf_df.agg(['mean', 'median', 'std', 'min', 'max'])
    logging.info(f"TF-IDF statistics:\n{stats}")
    
    # Plot boxplots of TF-IDF scores for each region
    plt.figure(figsize=(10, 6))
    tfidf_df.boxplot()
    plt.title("Distribution of TF-IDF Scores by Region")
    plt.ylabel("TF-IDF Score")
    plt.yscale('log')  # Use log scale for better visualization
    plt.tight_layout()
    plt.show()

def pairwise_wilcoxon(tfidf_df):
    """
    Perform pairwise Wilcoxon signed-rank tests between regions.
    Returns and logs a DataFrame of test results.
    """
    logging.info("Performing pairwise Wilcoxon signed-rank tests...")
    regions = tfidf_df.columns
    results = []
    
    for i in range(len(regions)):
        for j in range(i+1, len(regions)):
            region1, region2 = regions[i], regions[j]
            statistic, p_value = wilcoxon(tfidf_df[region1], tfidf_df[region2])
            results.append((region1, region2, statistic, p_value))
    
    results_df = pd.DataFrame(results, columns=['Region1', 'Region2', 'Statistic', 'p-value'])
    logging.info(f"Wilcoxon test results:\n{results_df}")
    return results_df

def top_distinctive_terms(tfidf_df, n=10):
    """
    Identify and log the top distinctive terms for each region.
    """
    logging.info(f"Identifying top {n} distinctive terms for each region...")
    
    for region in tfidf_df.columns:
        # Calculate the difference between this region's scores and the mean of other regions
        distinctiveness = tfidf_df[region] - tfidf_df.drop(columns=[region]).mean(axis=1)
        top_terms = distinctiveness.nlargest(n)
        
        logging.info(f"Top {n} distinctive terms for {region}:")
        for term, score in top_terms.items():
            logging.info(f"  {term}: {score:.4f}")

def calculate_expansion_weighted_scores(tfidf_df, coordinate_sets):
    """
    Calculate expansion-weighted scores for each term.
    Returns a Series of weighted scores.
    """
    expansion_factors = {region: data['scaling_factor'] for region, data in coordinate_sets.items()}
    weighted_scores = pd.Series(0, index=tfidf_df.index)
    
    for region in tfidf_df.columns:
        weighted_scores += tfidf_df[region] * expansion_factors[region]
    
    return weighted_scores

def generate_expansion_word_cloud(tfidf_df, coordinate_sets):
    """
    Generate a word cloud of terms associated with cortical expansion.
    Displays the word cloud and logs top terms.
    """
    logging.info("Generating expansion-associated word cloud...")
    
    weighted_scores = calculate_expansion_weighted_scores(tfidf_df, coordinate_sets)
    print(weighted_scores)
    
    # Normalize scores to be non-negative
    min_score = weighted_scores.min()
    if min_score < 0:
        weighted_scores = weighted_scores - min_score

     # Initialize a dictionary to count the number of regions each word appears in the top 30
    term_region_counts = {word: 0 for word in weighted_scores.index}
    
    # Iterate over each region to find the top 30 words and update the counts
    for region in tfidf_df.columns:
        top_30_words_in_region = tfidf_df[region].nlargest(50).index
        for word in top_30_words_in_region:
            if word in term_region_counts:
                term_region_counts[word] += 1
    
    print(term_region_counts)
    
    # Define a function to map the count to a color
    def term_color_func(word, *args, **kwargs):
        count = term_region_counts.get(word, 0)
        print(count)
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

def perform_chi_square_test(tfidf_df, percentile=50):
    """
    Perform a chi-square test on the TF-IDF data.
    Returns chi-square statistic and p-value if applicable.
    """
    logging.info(f"Performing chi-square test using {percentile}th percentile as threshold...")
    threshold = tfidf_df.quantile(percentile/100).quantile(percentile/100)
    contingency_table = (tfidf_df > threshold).astype(int)
    
    # Log information about the contingency table
    logging.info(f"Contingency table shape: {contingency_table.shape}")
    logging.info(f"Contingency table sum: {contingency_table.sum().sum()}")
    logging.info(f"Threshold used: {threshold}")
    
    # Check if any row or column sums to zero
    row_sums = contingency_table.sum(axis=1)
    col_sums = contingency_table.sum(axis=0)
    
    if (row_sums == 0).any():
        logging.warning(f"Number of zero rows: {(row_sums == 0).sum()}")
    if (col_sums == 0).any():
        logging.warning(f"Number of zero columns: {(col_sums == 0).sum()}")
    
    if (row_sums == 0).any() or (col_sums == 0).any():
        logging.warning("Contingency table has zero rows or columns. Chi-square test may not be appropriate.")
        return None, None
    
    try:
        chi2, p_value, dof, expected = chi2_contingency(contingency_table)
        logging.info(f"Chi-square statistic: {chi2}")
        logging.info(f"p-value: {p_value}")
        return chi2, p_value
    except ValueError as e:
        logging.error(f"Error in chi-square test: {e}")
        return None, None

def main():
    """
    Main function to orchestrate the entire analysis process.
    """
    metadata_df, coordinates_df, feature_df = load_data()
    coordinate_sets, new_data = define_coordinate_sets()
    print(new_data)
    np.savetxt("./rois/coords_w_expansion.csv", new_data, delimiter=",", fmt='%.3f')
    combined_term_data = process_terms(metadata_df, coordinates_df, feature_df, coordinate_sets)
    
    if combined_term_data.empty:
        logging.error("No terms extracted for any region. Exiting.")
        return

    tfidf_df = calculate_tfidf(combined_term_data)
    
    logging.info(f"TF-IDF DataFrame shape: {tfidf_df.shape}")
    logging.info(f"TF-IDF DataFrame columns: {tfidf_df.columns.tolist()}")
    
    generate_word_clouds(tfidf_df, coordinate_sets)
    generate_heatmap(tfidf_df)
    generate_bar_graph(tfidf_df)
    analyze_tfidf_distribution(tfidf_df)
    pairwise_wilcoxon(tfidf_df)
    top_distinctive_terms(tfidf_df)
    generate_expansion_word_cloud(tfidf_df, coordinate_sets)
    
    logging.info("Analysis complete.")

if __name__ == "__main__":
    main()