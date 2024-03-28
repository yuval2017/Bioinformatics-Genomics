import pandas as pd

path_5UTR = '../Data/5UTR_features2.csv'
path_CDS = '../Data/CDS_features2.csv'
path_3UTR = '../Data/3UTR_features2.csv'


class Data:
    def __init__(self):
        self.filter_column = ['ReadSequences', 'miRNA', 'Read_start_5', 'Read_end_5', 'mRNA', 'Read_start_3',
                              'Read_end_3',
                              'Binding_Region']

        save_csv_path = '../Data/Table_S2.csv'

        # Load data into a Pandas DataFrame
        df = pd.read_csv(save_csv_path)
        # filter columns
        self.df = df[self.filter_column]
