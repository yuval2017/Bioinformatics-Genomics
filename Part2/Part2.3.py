import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import f_oneway
from itertools import combinations
from Part2.Data import *

filter_columns2 = ["Energy", "BP", "GC_BP", "AU_BP", "GU_BP", "mismatches", "miRNA_b", 'mRNA_b', "Conservation Score",
                  "Normalized Accessibility", "MFE Value"]
filter_columns = ["Energy", "BP", "GC_BP", "AU_BP", "miRNA_b", "Conservation Score",
                  "Normalized Accessibility"]

def plot_graphs(region_pds):
    region_pairs = region_pds.items()
    # create tabula
    fig, axs = plt.subplots(len(filter_columns), 3, figsize=(12, 18))
    for i, feature in enumerate(filter_columns):
        for j, (region, pd_region) in enumerate(region_pairs):
            pd_region = pd_region.dropna(subset=[feature])

            pd_region[feature] = pd_region[feature].astype(float)

            pd_region = pd_region.sort_values(by=feature)
            axs[i, j].hist(pd_region[feature], bins=30, color='red', edgecolor='black', alpha=0.7)
            axs[i, j].set_xlabel('')
            axs[i, j].set_ylabel('')
            axs[i, j].set_title(f'{feature} - {region}')

    plt.tight_layout(rect=[0, 0, 1, 1], h_pad=1.5, w_pad=2)
    plt.show()


def anova_test(df, regions, features):
    # Group your data by 'region'
    grouped_data = [df[df['region'] == region][features] for region in regions]

    # Perform ANOVA test
    f_statistic, p_value = f_oneway(*grouped_data)

    # Output results
    print("ANOVA Test Results for:", *regions)
    print("F-statistic:", f_statistic)
    print("P-value:", p_value)

    # Plot boxplots for visual comparison

    plt.boxplot(grouped_data, labels=regions)
    plt.title("Distribution of Feature Across Regions")
    plt.xlabel("Region")
    plt.ylabel("Feature Value")
    plt.show()
def find_subgroups(strings):
    subgroups = []
    for r in range(2, len(strings) + 1):  # Iterate over subgroup sizes from 2 to length of array
        for combo in combinations(strings, r):
            subgroups.append(combo)
    return subgroups


if __name__ == "__main__":
    find_subgroups(['3UTR', '5UTR', 'CDS'])
    pd_5UTR = pd.read_csv(path_5UTR)
    pd_CDS = pd.read_csv(path_CDS)
    pd_3UTR = pd.read_csv(path_3UTR)


    # ass 2.3.1
    region_pds = {"5UTR": pd_5UTR, "CDS": pd_CDS, "3UTR": pd_3UTR}
    plot_graphs(region_pds)

    # ass 2.3.2
    pd_5UTR['region'] = '5UTR'
    pd_CDS['region'] = 'CDS'
    pd_3UTR['region'] = '3UTR'
    combined_data = pd.concat([pd_5UTR, pd_CDS, pd_3UTR], ignore_index=True)
    feature = 'AU_BP'
    anova_test(combined_data, ['3UTR', '5UTR', 'CDS'], feature)
    anova_test(combined_data, ['3UTR', '5UTR'], feature)
    anova_test(combined_data, ['3UTR', 'CDS'], feature)
    anova_test(combined_data, ['5UTR', 'CDS'], feature)
