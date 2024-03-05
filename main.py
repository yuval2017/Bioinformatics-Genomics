import random


# Function to generate a random DNA sequence of given length
def generate_sequence(length):
    nucleotides = ['A', 'C', 'T', 'G']
    sequence = [random.choice(nucleotides) for _ in range(length)]
    return ''.join(sequence)


# Generate random sequences
def generate_sequences(length):
    # Generate the first sequence
    s1 = generate_sequence(length)

    # Generate the second sequence ensuring that each nucleotide appears at least once
    s2 = generate_sequence(length - 1)  # Start with one less nucleotide
    missing_nucleotide = random.choice([nt for nt in 'ACTG' if nt not in s2])
    s2 += missing_nucleotide  # Add the missing nucleotide to the second sequence

    return s1, s2


def find_global_alignments_with_scores(seq1, seq2, match=1, mismatch=-1, indel=-1, max_alignments=2):
    n, m = len(seq1), len(seq2)
    dp = [[0 for _ in range(m + 1)] for _ in range(n + 1)]

    # Initialize first row and column for global alignment
    for i in range(1, n + 1):
        dp[i][0] = i * indel
    for j in range(1, m + 1):
        dp[0][j] = j * indel

    # Fill DP table
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i - 1] == seq2[j - 1]:
                dp[i][j] = max(dp[i - 1][j - 1] + match, dp[i - 1][j] + indel, dp[i][j - 1] + indel)
            else:
                dp[i][j] = max(dp[i - 1][j - 1] + mismatch, dp[i - 1][j] + indel, dp[i][j - 1] + indel)

    def traceback(i, j, align1="", align2="", trace=[], paths=[]):
        if len(paths) >= max_alignments:
            return
        if i > 0 and j > 0:
            if seq1[i - 1] == seq2[j - 1]:
                if dp[i][j] == dp[i - 1][j - 1] + match:
                    traceback(i - 1, j - 1, seq1[i - 1] + align1, seq2[j - 1] + align2, [(i, j)] + trace, paths)
            else:
                if dp[i][j] == dp[i - 1][j - 1] + mismatch:
                    traceback(i - 1, j - 1, seq1[i - 1] + align1, seq2[j - 1] + align2, [(i, j)] + trace, paths)
        if i > 0 and dp[i][j] == dp[i - 1][j] + indel:
            traceback(i - 1, j, seq1[i - 1] + align1, "-" + align2, [(i, j)] + trace, paths)
        if j > 0 and dp[i][j] == dp[i][j - 1] + indel:
            traceback(i, j - 1, "-" + align1, seq2[j - 1] + align2, [(i, j)] + trace, paths)
        if i == 0 and j == 0:
            paths.append((align1, align2, trace[::-1]))  # Add the alignment and its trace reversed for correct order

    paths = []
    traceback(n, m, "", "", [], paths)

    alignments_with_traces_and_scores = []
    for align1, align2, trace in paths:
        score = dp[n][m]  # The score of global alignment is in the bottom-right cell
        alignments_with_traces_and_scores.append({"alignment": (align1, align2), "trace": trace, "score": score})

    return {"alignments_with_traces_and_scores": alignments_with_traces_and_scores, "dp_table": dp,
            "optimal_score": dp[n][m]}


# Re-import NumPy and define the function again
import numpy as np


def ends_space_free_alignment_with_ends(seq1, seq2, match_score=1, mismatch_score=-1, gap_score=-1):
    rows, cols = len(seq1) + 1, len(seq2) + 1
    score_matrix = np.zeros((rows, cols), dtype=int)

    # Fill the scoring matrix
    for i in range(1, rows):
        for j in range(1, cols):
            if seq1[i-1] == seq2[j-1]:
                match = score_matrix[i-1, j-1] + match_score
            else:
                match = score_matrix[i-1, j-1] + mismatch_score
            delete = score_matrix[i-1, j] + gap_score
            insert = score_matrix[i, j-1] + gap_score
            score_matrix[i, j] = max(match, delete, insert, 0)

    # Traceback to find the alignment, including ends
    def traceback(start_pos):
        i, j = start_pos
        align1, align2 = "", ""
        path = [(i, j)]
        while i > 0 or j > 0:
            if i > 0 and (j == 0 or score_matrix[i, j] == score_matrix[i-1, j] + gap_score):
                align1 = seq1[i-1] + align1
                align2 = "-" + align2
                i -= 1
            elif j > 0 and (i == 0 or score_matrix[i, j] == score_matrix[i, j-1] + gap_score):
                align1 = "-" + align1
                align2 = seq2[j-1] + align2
                j -= 1
            else:
                align1 = seq1[i-1] + align1
                align2 = seq2[j-1] + align2
                i -= 1
                j -= 1
            path.insert(0, (i, j))

        # Add trailing '-' to align ends if necessary
        while len(align1) < len(align2):
            align1 += "-"
        while len(align2) < len(align1):
            align2 += "-"

        return align1, align2, path, score_matrix[start_pos[0], start_pos[1]]

    # Find the starting positions for the highest scoring alignments
    best_scores = np.argwhere(score_matrix == np.max(score_matrix))
    alignments = []
    for start_pos in best_scores[:2]:  # Up to two alignments
        alignments.append(traceback((start_pos[0], start_pos[1])))

    return alignments, score_matrix

# Example use


# Example usage

# Example usage
# Test the function with example sequences


if __name__ == "__main__":
    # Generate sequences of length 10
    #    s1, s2 = generate_sequences(10)

    # Print the generated sequences
    # print("Sequence 1 (S1):", s1)
    # print("Sequence 2 (S2):", s2)

    seq1 = "CACTGTAC"
    seq2 = "GACACTTG"
    result = find_global_alignments_with_scores(seq1, seq2, 2)
    for item in result["alignments_with_traces_and_scores"]:
        print("Alignment:", item["alignment"])
        print("Trace:", item["trace"])
        print("Score:", item["score"])
        print()

    print('--------------------')
    alignments_with_ends, score_matrix_with_ends = ends_space_free_alignment_with_ends("CACTGTAC", "GACACTTG", 2)
    print(alignments_with_ends)
    print(score_matrix_with_ends)

    # Example usage

    # dp_table, aligned_seq1, aligned_seq2 = end_space_free_alignment(seq1, seq2)
