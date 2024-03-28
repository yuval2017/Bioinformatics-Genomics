import pandas as pd

features = ["miRNA", "mRNA", "miRNA fragment (5\'-3\')", "mRNA fragment (3\'-5\')", "RNAduplex output", "miRNA interacting region(5'-3')", "mRNA interacting region(3'-5')", "Duplex", "Energy", "BP", "GC_BP",
            "AU_BP", "miRNA_b", "Conservation Score", "Normalized Accessibility"]

if __name__ == "__main__":
    save_csv_path = '../Data/Table_S2.csv'

    df_5UTR = pd.read_csv('../Data/5UTR_features2.csv')
    df_CDS = pd.read_csv('../Data/CDS_features2.csv')
    df_3UTR = pd.read_csv('../Data/3UTR_features2.csv')

    df_5UTR = df_5UTR[features]
    df_CDS = df_CDS[features]
    df_3UTR = df_3UTR[features]

    df_5UTR.to_csv('../Data/5UTR_features_sub.csv', index=False)
    df_CDS.to_csv('../Data/CDS_features_sub.csv', index=False)
    df_3UTR.to_csv('../Data/3UTR_features_sub.csv', index=False)
