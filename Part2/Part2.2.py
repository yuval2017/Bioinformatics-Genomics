import re
import pandas as pd
import subprocess
from Bio import pairwise2
from Bio.Align import substitution_matrices
import RNA
def extract_miRNA_fragment_5(row):
    read_sequence = row['ReadSequences']
    start_index = row['Read_start_5']
    end_index = row['Read_end_5']
    return read_sequence[start_index - 1:end_index]


def extract_mRNA_fragment_3(row):
    read_sequence = row['ReadSequences']
    start_index = row['Read_start_3']
    end_index = row['Read_end_3']
    return read_sequence[start_index - 1:end_index]


def RNA_duplex(miRNA, mRNA):
    try:
        # Run RNAduplex
        command = subprocess.Popen(['RNAduplex'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                   universal_newlines=True)
        input_seq = f"{miRNA}\n{mRNA}\n"
        output, error = command.communicate(input=input_seq)

        value = output.split('\n')[0]
        val1, val2 = value.split(':')
        pattern1 = r"\s*([\.(\)&]+)\s+(\d+\s*\,\s*\d+)\s*"
        pattern2 = r"\s*(\d+\s*\,\s*[0-9]+)\s+\(\s*(-?\d+\.\d+)\s*\)\s*"
        # Match the pattern in the string
        matches1 = re.match(pattern1, val1)
        structure, pair_miRNA = matches1.groups()
        matches2 = re.match(pattern2, val2)
        pair_mRNA, energy = matches2.groups()
        if len(matches1.groups()) > 2 or len(matches2.groups()) > 2:
            print("Regular exp problem!")

        start_miRNA, end_miRNA = map(lambda x: int(x.strip()), pair_miRNA.split(','))
        start_mRNA, end_mRNA = map(lambda x: int(x.strip()), pair_mRNA.split(','))

        # Return parsed data
        return value, structure, energy, start_miRNA, end_miRNA, start_mRNA, end_mRNA
    except Exception as e:
        # Handle exceptions
        print("Error:", e)
        return None


def duplex_string(structure, miRNA, mRNA):
    structure_miRNA, structure_mRNA = structure.split('&')
    if len(miRNA) != len(structure_miRNA):
        print("here1")
    if len(mRNA) != len(structure_mRNA):
        print("here2")
    stack_miRNA = list(zip(list(miRNA), list(structure_miRNA)))[::-1]
    stack_mRNA = list(zip(list(mRNA), list(structure_mRNA)))[::-1]

    b_miRNA = b_mRNA = mismatch = matches = bp_gc = bg_au = bg_gu = 0
    top_line = middle_line1 = middle_line2 = bottom_line = ''
    while stack_miRNA and stack_mRNA:
        c_top = c_middle1 = c_middle2 = c_bottom = ' '

        if stack_miRNA[-1][1] == '.' and stack_mRNA[-1][1] == '.':
            c_top, _ = stack_miRNA.pop()
            c_bottom, _ = stack_mRNA.pop()
            mismatch += 1

        elif stack_miRNA[-1][1] == '.' and stack_mRNA[-1][1] == ')':
            c_top, _ = stack_miRNA.pop()
            b_miRNA += 1
        elif stack_miRNA[-1][1] == '(' and stack_mRNA[-1][1] == '.':
            c_bottom, _ = stack_mRNA.pop()
            b_mRNA += 1
        else:
            matches += 1
            c_middle1, _ = stack_miRNA.pop()
            c_middle2, _ = stack_mRNA.pop()
            if c_middle1 + c_middle2 in ['GC', 'CG']:
                bp_gc += 1
            if c_middle1 + c_middle2 in ['AT', 'TA']:
                bg_au += 1
            if c_middle1 + c_middle2 in ['GT', 'TG']:
                bg_gu += 1

        top_line += c_top
        middle_line1 += c_middle1
        middle_line2 += c_middle2
        bottom_line += c_bottom

    while stack_miRNA:
        if stack_miRNA[-1][1] != '.':
            print("something wrong1")

        b_miRNA += 1
        top_line += ' '
        middle_line1 += stack_miRNA.pop()[0]
        middle_line2 += ' '
        bottom_line += ' '

    while stack_mRNA:
        if stack_mRNA[-1][1] != '.':
            print("something wrong2")
        b_mRNA += 1
        top_line += ' '
        middle_line1 += ' '
        middle_line2 += stack_mRNA.pop()[0]
        bottom_line += ' '

    duplex = top_line + '\n' + middle_line1 + '\n' + middle_line2 + '\n' + bottom_line
    return duplex, matches, bp_gc, bg_au, bg_gu, mismatch, b_miRNA, b_mRNA


# Define a function to compute the RNA duplex
def compute_RNA_duplex(row):
    miRNA = row['miRNA fragment (5\'-3\')']
    mRNA = row['mRNA fragment (3\'-5\')']
    seq = row["ReadSequences"]
    out, structure, energy, start_miRNA, end_miRNA, start_mRNA, end_mRNA = RNA_duplex(miRNA, mRNA)

    if structure == '.&.':
        return out, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"

    miRNA_interacting_region = miRNA[start_miRNA - 1: end_miRNA]
    mRNA_interacting_region = mRNA[start_mRNA - 1: end_mRNA][::-1]
    duplex, matches, bp_gc, bg_au, bg_gu, mismatch, b_miRNA, b_mRNA = duplex_string(structure, miRNA_interacting_region,
                                                              mRNA_interacting_region)
    conservation_score = sequence_conservation_score(miRNA, mRNA)
    normalized_accessibility, mfe = calculate_accessibility(seq)
    return out, miRNA_interacting_region, mRNA_interacting_region, duplex, energy, matches, bp_gc, bg_au, bg_gu, mismatch, b_miRNA, b_mRNA , conservation_score, normalized_accessibility, mfe

def sequence_conservation_score(miRNA_sequence, mRNA_sequence):
    # Perform global pairwise alignment using a DNA substitution matrix
    alignments = pairwise2.align.globalds(miRNA_sequence, mRNA_sequence, substitution_matrices.load("BLOSUM62"), -10, -0.5)

    # Calculate sequence identity as the number of matches divided by the alignment length
    alignment = alignments[0]
    alignment_length = len(alignment[0])
    num_matches = sum(1 for i in range(alignment_length) if alignment[0][i] == alignment[1][i])
    sequence_identity = num_matches / alignment_length
    return sequence_identity



def exec_and_save_df(df, csv_path):
    filter_columns = ["miRNA", "mRNA", "miRNA fragment (5\'-3\')", "mRNA fragment (3\'-5\')", "RNAduplex output",
                      "miRNA interacting region(5'-3')", "mRNA interacting region(3'-5')", "Duplex", "Energy", "BP",
                      "GC_BP","AU_BG", "GU_BP", "mismatches", "miRNA_b", 'mRNA_b', "Conservation Score", "Normalized Accessibility", "MFE Value"]
    # Apply the function to each row and create a new column
    df['miRNA fragment (5\'-3\')'] = df.apply(extract_miRNA_fragment_5, axis=1)
    df['mRNA fragment (3\'-5\')'] = df.apply(extract_mRNA_fragment_3, axis=1)

    # Apply the function to each row and create a new column
    df['RNAduplex output'] = df.apply(compute_RNA_duplex, axis=1)
    df['RNAduplex output'], df["miRNA interacting region(5'-3')"], df["mRNA interacting region(3'-5')"], df['Duplex'], \
    df['Energy'], df['BP'], df['GC_BP'], df['AU_BG'],  df['GU_BP'], df['mismatches'], df["miRNA_b"], df['mRNA_b'], df["Conservation Score"], df["Normalized Accessibility"], df["MFE Value"] = zip(*df.apply(compute_RNA_duplex, axis=1))
    df = df[filter_columns]
    # Save the DataFrame to a CSV file
    df.to_csv(csv_path, index=False)

def calculate_accessibility(rna_sequence):
    # Predict the minimum free energy (MFE) secondary structure
    (structure, mfe) = RNA.fold(rna_sequence)

    # Calculate the accessibility based on the base pairs in the secondary structure
    accessibility = 0
    for i in range(len(structure)):
        if structure[i] == '.':
            accessibility += 1

    # Normalize the accessibility score by the length of the sequence
    normalized_accessibility = accessibility / len(rna_sequence)

    return normalized_accessibility, mfe

from collections import Counter

def calculate_kmer_features(rna_sequence, k):
    # Initialize a Counter to store the frequency of each k-mer
    kmer_counter = Counter()

    # Iterate over the RNA sequence to extract k-mers
    for i in range(len(rna_sequence) - k + 1):
        # Extract the k-mer of length k starting at position i
        kmer = rna_sequence[i:i + k]
        # Increment the count for this k-mer
        kmer_counter[kmer] += 1

    # Normalize the counts to obtain frequencies
    total_kmers = sum(kmer_counter.values())
    kmer_features = {kmer: count / total_kmers for kmer, count in kmer_counter.items()}

    return kmer_features

def most_frequent_kmer(df, k):
    # Initialize a Counter to store the occurrence counts of k-mers
    kmer_counts = Counter()

    # Iterate over each sequence in the DataFrame
    for sequence in df['ReadSequences']:
        # Calculate the occurrence of k-mers for the current sequence
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            kmer_counts[kmer] += 1

    # Find the most frequent k-mer
    most_common_kmer = kmer_counts.most_common(1)

    return most_common_kmer

if __name__ == "__main__":
    # Example usage


    #
    # rna_sequence = "AUUGCCAGUUCGAUUCGGAAUUCGU"
    #
    # k = 3  # Length of k-mer
    # #kmer_features = calculate_kmer_features(rna_sequence, k)
    # accessibility = calculate_accessibility(rna_sequence)
    # seq_miRNA = "AAGCTGCCAGTTGAAGAACTGT"
    # seq_mRNA = "TTACTTCATGGCAGCTATCCCACAG"
    # ans1 = sequence_conservation_score(seq_miRNA, seq_mRNA)
    # print(ans1)
    # #ans2 = calculate_conservation_score([seq_miRNA, seq_mRNA])
    # # Your string
    # value = ".((((((((..(((((.&.))))))))))))).   1,17  :   3,17  (-19.70)"
    # val1, val2 = value.split(':')
    # pattern1 = r"\s*([\.(\)&]+)\s+(\d+\s*\,\s*\d+)\s*"
    # pattern2 = r"\s*(\d+\s*\,\s*[0-9]+)\s+\(\s*(-?\d+\.\d+)\s*\)\s*"
    # # Match the pattern in the string
    # matches1 = re.match(pattern1, val1)
    # grp1 = matches1.groups()
    # matches2 = re.match(pattern2, val2)
    # grp2 = matches2.groups()
    # duplex_string('.((((((((..(((((.&.))))))))))))).', 'AAGCTGCCAGTTGAAGA', 'ATCGACGGTACTTCA')
    # save_csv_path = '../Data/Table_S2.csv'
    # example = 'AAGCTGCCAGTTGAAGAACTGTTTACTTCATGGCAGCTATCCCACA'
    # print(example[22:47])
    # Load data into a Pandas DataFrame

    save_csv_path = '../Data/Table_S2.csv'
    df = pd.read_csv(save_csv_path)
    df_5UTR = df[df['Binding_Region'] == "5'UTR"].copy()
    df_CDS = df[df['Binding_Region'] == 'CDS'].copy()
    df_3UTR = df[df['Binding_Region'] == "3'UTR"].copy()

    exec_and_save_df(df_5UTR, '../Data/5UTR_features2.csv')
    exec_and_save_df(df_CDS, '../Data/CDS_features2.csv')
    exec_and_save_df(df_3UTR, '../Data/3UTR_features2.csv')
