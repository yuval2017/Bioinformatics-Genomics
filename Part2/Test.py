import subprocess
def analyze_structure(miRNA_fragment, mRNA_fragment, structure):
    # Split the structure at the '&' to analyze miRNA and mRNA separately
    miRNA_structure, mRNA_structure = structure.split('&')

    # Initialize counters for matched bases in miRNA and mRNA
    miRNA_matches = [0] * len(miRNA_fragment)
    mRNA_matches = [0] * len(mRNA_fragment)

    # Initialize variables to hold the pairing information
    stack = []
    pair_map = {}

    # Analyze the miRNA structure for base pairing
    for i, char in enumerate(miRNA_structure):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                opening = stack.pop()
                closing = i
                pair_map[opening] = closing
                miRNA_matches[opening] = 1
                miRNA_matches[closing] = 1

    # Clear the stack for mRNA analysis
    stack.clear()

    # Analyze the mRNA structure for base pairing
    for i, char in enumerate(mRNA_structure):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                opening = stack.pop()
                closing = i
                # Adjust index for combined structure
                adjusted_opening = len(miRNA_fragment) + 1 + opening
                adjusted_closing = len(miRNA_fragment) + 1 + closing
                pair_map[adjusted_opening] = adjusted_closing
                mRNA_matches[opening] = 1
                mRNA_matches[closing] = 1

    # Print base pairing information
    print("miRNA base pairings:", miRNA_matches)
    print("mRNA base pairings:", mRNA_matches)
    print("Pair map:", pair_map)
def predict_RNA_duplex(miRNA, mRNA):
    command = subprocess.Popen(['RNAduplex'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               universal_newlines=True)
    input_seq = f"{miRNA}\n{mRNA}\n"
    output, error = command.communicate(input=input_seq)
if __name__ == "__main__":
    # Example usage with dummy data
    miRNA_fragment = "AAGCTGCCAGTTGAAGAACTGT"
    mRNA_fragment = "TTACTTCATGGCAGCTATCCCACAG"
    predict_RNA_duplex(miRNA_fragment, mRNA_fragment)
    predicted_structure = ".((((((((..(((((.&.)))))))))))))."

    analyze_structure(miRNA_fragment, mRNA_fragment, predicted_structure)