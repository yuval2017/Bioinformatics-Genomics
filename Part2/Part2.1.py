import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3


if __name__ == "__main__":
    filter_column = ['ReadSequences', 'miRNA', 'Read_start_5', 'Read_end_5', 'mRNA', 'Read_start_3', 'Read_end_3',
                     'Binding_Region']

    # Define file paths

    save_csv_path = '../Data/Table_S2.csv'

    # Load data into a Pandas DataFrame
    df = pd.read_csv(save_csv_path)
    # filter columns
    df = df[filter_column]

    # Define the values you want to keep
    values_to_keep = ["5'UTR", "CDS", "3'UTR"]
    for region in values_to_keep:
        print(f"Total number of interactions per miRNA and mRNA in region {region}:",
              len(df[df["Binding_Region"] == region]))
    # Filter the DataFrame to keep rows with specified values
    filtered_df = df[df['Binding_Region'].isin(values_to_keep)]

    counts_miRNA = filtered_df.groupby(['miRNA', 'Binding_Region']).size().reset_index(name='count')
    counts_mRNA = filtered_df.groupby(['mRNA', 'Binding_Region']).size().reset_index(name='count')

    # Pivot the table to create a new DataFrame with miRNA as rows, Binding_Region as columns, and counts as values
    counts_pivot_miRNA = counts_miRNA.pivot_table(index='miRNA', columns='Binding_Region', values='count', fill_value=0)
    counts_pivot_mRNA = counts_mRNA.pivot_table(index='mRNA', columns='Binding_Region', values='count', fill_value=0)

    top_20_miRNAs_indexes = counts_pivot_miRNA.sum(axis=1).nlargest(20).index
    top_20_mRNAs_indexes = counts_pivot_mRNA.sum(axis=1).nlargest(20).index

    top_20_counts_miRNA = counts_pivot_miRNA.loc[top_20_miRNAs_indexes]
    top_20_counts_mRNA = counts_pivot_mRNA.loc[top_20_mRNAs_indexes]

    bar_colors = ['blue', 'green', 'red']
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))

    #top_20_counts_miRNA.plot(kind='bar', ax=axes[0], stacked=True, color=bar_colors)
    #top_20_counts_mRNA.plot(kind='bar', ax=axes[1], stacked=True, color=bar_colors)

    # # Set titles and axis labels
    # axes[0].set_title('Top 20 miRNA Interactions by Binding Region')
    # axes[0].set_xlabel('miRNA')
    # axes[0].set_ylabel('Total Number Of Interactions')
    #
    # axes[1].set_title('Top 20 mRNA Interactions by Binding Region')
    # axes[1].set_xlabel('mRNA')
    # axes[1].set_ylabel('Total Number Of Interactions')
    #
    # plt.tight_layout()  # Adjust layout to prevent overlap of labels
    # plt.show()
    # plt.subplots(figsize=(12, 10))
    # top_miRNA_5UTR = set(counts_pivot_miRNA.nlargest(20, "5'UTR").index)
    # top_miRNA_3UTR = set(counts_pivot_miRNA.nlargest(20, "3'UTR").index)
    # top_miRNA_CDS = set(counts_pivot_miRNA.nlargest(20, "CDS").index)
    # venn3([top_miRNA_5UTR, top_miRNA_3UTR, top_miRNA_CDS], values_to_keep)
    # plt.title("Overlap between top miRNAs at different regions")
    # plt.show()

