import re

import pandas as pd
from ViennaRNA import RNA
import subprocess


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
        # Run RNAduplex subprocess
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
        if len(matches1.groups()) > 2 or  len(matches2.groups()) > 2:
            print("Regular exp problem!")

        start_miRNA, end_miRNA = map(lambda x: int(x.strip()), pair_miRNA.split(','))
        start_mRNA, end_mRNA = map(lambda x: int(x.strip()), pair_mRNA.split(','))

        # Return parsed data
        return value, structure, energy, start_miRNA, end_miRNA, start_mRNA, end_mRNA
    except Exception as e:
        # Handle exceptions
        print("Error:", e)
        return None


# Define a function to compute the RNA duplex
def compute_RNA_duplex(row):
    miRNA = row['miRNA fragment (5\'-3\')']
    mRNA = row['mRNA fragment (3\'-5\')']
    out, structure, energy, start_miRNA, end_miRNA, start_mRNA, end_mRNA = RNA_duplex(miRNA, mRNA)

    if structure == '.&.':
        return out, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"

    miRNA_interacting_region = miRNA[start_miRNA - 1: end_miRNA]
    mRNA_interacting_region = mRNA[start_mRNA - 1: end_mRNA][::-1]
    duplex, matches, bp_gc, mismatch, b_miRNA = duplex_string(structure, miRNA_interacting_region,
                                                              mRNA_interacting_region)
    return out, miRNA_interacting_region, mRNA_interacting_region, duplex, energy, matches, bp_gc, mismatch, b_miRNA


def duplex_string(structure, miRNA, mRNA):
    structure_miRNA, structure_mRNA = structure.split('&')
    if len(miRNA) != len(structure_miRNA):
        print("here1")
    if len(mRNA) != len(structure_mRNA):
        print("here2")
    stack_miRNA = list(zip(list(miRNA), list(structure_miRNA)))[::-1]
    stack_mRNA = list(zip(list(mRNA), list(structure_mRNA)))[::-1]

    b_miRNA = b_mRNA = mismatch = matches = bp_gc = 0
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

    return duplex, matches, bp_gc, mismatch, b_miRNA


def exec_and_save_df(df, csv_path):
    filter_columns = ["miRNA", "mRNA", "miRNA fragment (5\'-3\')", "mRNA fragment (3\'-5\')", "RNAduplex output",
                      "miRNA interacting region(5'-3')", "mRNA interacting region(3'-5')", "Duplex", "Energy", "BP",
                      "GC_BP", "mismatches", "miRNA_b"]
    # Apply the function to each row and create a new column
    df['miRNA fragment (5\'-3\')'] = df.apply(extract_miRNA_fragment_5, axis=1)
    df['mRNA fragment (3\'-5\')'] = df.apply(extract_mRNA_fragment_3, axis=1)

    # Apply the function to each row and create a new column
    df['RNAduplex output'] = df.apply(compute_RNA_duplex, axis=1)
    df['RNAduplex output'], df["miRNA interacting region(5'-3')"], df["mRNA interacting region(3'-5')"], df['Duplex'], \
    df['Energy'], df['BP'], df['GC_BP'], df['mismatches'], df["miRNA_b"] = zip(*df.apply(compute_RNA_duplex, axis=1))
    df = df[filter_columns]
    # Save the DataFrame to a CSV file
    df.to_csv(csv_path, index=False)


if __name__ == "__main__":
    miRNA = 20 * 'A'
    mRNA = 20 * 'C'
    duplex1 = RNA.duplexfold(miRNA, mRNA)
    structure = duplex1.structure
    energy = duplex1.energy
    duplex2 = RNA.duplexfold(mRNA, miRNA)
    # Your string
    value = ".((((((((..(((((.&.))))))))))))).   1,17  :   3,17  (-19.70)"
    val1, val2 = value.split(':')
    pattern1 = r"\s*([\.(\)&]+)\s+(\d+\s*\,\s*\d+)\s*"
    pattern2 = r"\s*(\d+\s*\,\s*[0-9]+)\s+\(\s*(-?\d+\.\d+)\s*\)\s*"
    # Match the pattern in the string
    matches1 = re.match(pattern1, val1)
    grp1 = matches1.groups()
    matches2 = re.match(pattern2, val2)
    grp2 = matches2.groups()
    duplex_string('.((((((((..(((((.&.))))))))))))).', 'AAGCTGCCAGTTGAAGA', 'ATCGACGGTACTTCA')
    save_csv_path = '../Data/Table_S2.csv'
    example = 'AAGCTGCCAGTTGAAGAACTGTTTACTTCATGGCAGCTATCCCACA'
    print(example[22:47])
    # Load data into a Pandas DataFrame
    df = pd.read_csv(save_csv_path)

    df_5UTR = df[df['Binding_Region'] == "5'UTR"].copy()
    df_CDS = df[df['Binding_Region'] == 'CDS'].copy()
    df_3UTR = df[df['Binding_Region'] == "3'UTR"].copy()
    exec_and_save_df(df_5UTR, './Data/5UTR_features')
    exec_and_save_df(df_CDS, './Data/CDS_features')
    exec_and_save_df(df_3UTR, './Data/3UTR_features')
