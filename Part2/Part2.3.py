import pandas as pd
from matplotlib import pyplot as plt

path_5UTR = '../Data/5UTR_features.csv'
path_CDS = '../Data/CDS_features.csv'
path_3UTR = '../Data/3UTR_features.csv'
filter_columns = ["Energy", "BP", "GC_BP", "mismatches", "miRNA_b", "Conservation Score", "Normalized Accessibility", "MFE Value"]


if __name__ == "__main__":
    pd_5UTR = pd.read_csv(path_5UTR)
    pd_CDS = pd.read_csv(path_CDS)
    pd_3UTR = pd.read_csv(path_3UTR)
    region_pds = {"5UTR": pd_5UTR, "CDS": pd_CDS, "3UTR": pd_3UTR}
    region_pairs = region_pds.items()
    # create tabula
    fig, axs = plt.subplots(len(filter_columns), 3, figsize=(12, 18))
    for i, feature in enumerate(filter_columns):
        for j, (region, pd_region) in enumerate(region_pairs):
            pd_region = pd_region.dropna(subset=[feature])


            pd_region[feature] = pd_region[feature].astype(float)

            pd_region = pd_region.sort_values(by=feature)
            axs[i, j].hist(pd_region[feature], bins=20, color='red', edgecolor='black', alpha=0.7)
            axs[i, j].set_xlabel('')
            axs[i, j].set_ylabel('')
            axs[i, j].set_title(f'{feature} - {region}')

    plt.tight_layout(rect=[0, 0, 1, 0.95], h_pad=1.5, w_pad=2)
    plt.show()