import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

# Use the venn2 function

if __name__ == "__main__":
    filter_column = ['ReadSequences', 'miRNA', 'Read_start_5', 'Read_end_5', 'mRNA', 'Read_start_3', 'Read_end_3',
                     'Binding_Region']

    # Define file paths
    excel_file_path = '../Data/Table_S2.xlsx'
    save_csv_path = '../Data/Table_S2.csv'

    # Load data into a Pandas DataFrame
    df = pd.read_csv(save_csv_path)
    # filter columns
    df = df[filter_column]

    binding_regions = ['5’UTR', 'CDS', '3’UTR']
    filtered_df = df[df['Binding_Region'].isin(binding_regions)]

    # Compute the total number of interactions per miRNA and mRNA per region
    interactions_per_region = filtered_df.groupby(['miRNA', 'mRNA', 'Binding_Region']).size().reset_index(
        name='Total_Interactions')

    # Sum up the total interactions for each selected miRNA and mRNA
    total_interactions_per_miRNA = interactions_per_region.groupby('miRNA')['Total_Interactions'].sum()
    total_interactions_per_mRNA = interactions_per_region.groupby('mRNA')['Total_Interactions'].sum()


    top_20_interactions_per_miRNA = total_interactions_per_miRNA.nlargest(20)
    top_20_interactions_per_mRNA = total_interactions_per_mRNA.nlargest(20)
    print("Total Interactions per Top 20 miRNAs:")
    print(top_20_interactions_per_miRNA)
    print("\nTotal Interactions per Top 20 mRNAs:")
    print(top_20_interactions_per_mRNA)

    # # Plot total interactions for top 20 miRNAs
    # plt.figure(figsize=(12, 6))
    # top_20_interactions_per_miRNA.plot(kind='bar', color='blue')
    # plt.title('Total Interactions per Top 20 miRNAs')
    # plt.xlabel('miRNA')
    # plt.ylabel('Total Interactions')
    # plt.xticks(rotation=45, ha='right')
    # plt.tight_layout()
    # plt.show()
    #
    # # Plot total interactions for top 20 mRNAs
    # plt.figure(figsize=(12, 6))
    # top_20_interactions_per_mRNA.plot(kind='bar', color='green')
    # plt.title('Total Interactions per Top 20 mRNAs')
    # plt.xlabel('mRNA')
    # plt.ylabel('Total Interactions')
    # plt.xticks(rotation=45, ha='right')
    # plt.tight_layout()
    # plt.show()

    # Group by miRNA and Binding_Region and count the interactions
    interactions_per_miRNA_region = filtered_df.groupby(['Binding_Region', 'miRNA']).size()

    # Print the number of interactions for each miRNA and region
    print("Number of Interactions per miRNA and Region:")
    print(interactions_per_miRNA_region)
