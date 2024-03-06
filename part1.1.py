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

    def traceback(i, j, align1="", align2="", trace=None, paths=None):
        if paths is None:
            paths = []
        if trace is None:
            trace = []
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
    format = []
    for align1, align2, trace in paths:
        format.append((trace, align1, align2))

    return dp, dp[n][m], format


import numpy as np


def end_space_free_alignment(seq1, seq2, match_score=1, mismatch_score=-1, gap_score=-1):
    # Initialize the scoring matrix
    n, m = len(seq1), len(seq2)
    score_matrix = np.zeros((n + 1, m + 1), dtype=int)

    # Initialize first row and column based on gap penalty
    # For end-space free alignment, the first row and column are initialized to 0
    # as we allow free gaps at the beginning and end of the alignment

    # Fill in the scoring matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = score_matrix[i - 1, j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
            delete = score_matrix[i - 1, j] + gap_score
            insert = score_matrix[i, j - 1] + gap_score
            score_matrix[i, j] = max(match, delete, insert, 0)  # Allow dropping to 0 for local alignment

    # Traceback from the highest scoring cell
    # For end-space free alignment, we find the highest score in the last row or column
    max_score = np.max(score_matrix)
    max_pos = np.where(score_matrix == max_score)
    max_i, max_j = max_pos[0][0], max_pos[1][0]  # Take the first occurrence of the max score

    paths = [(max_i, max_j, [], [], [])]
    align_with_paths = []
    while paths:
        i, j, align1, align2, curr_path = paths.pop()
        curr_path = [(i, j)] + curr_path
        if i > 0 and j > 0 and score_matrix[i, j] > 0:
            score_current = score_matrix[i, j]
            score_diagonal = score_matrix[i - 1, j - 1]
            score_up = score_matrix[i - 1, j]
            score_left = score_matrix[i, j - 1]

            if score_current == score_diagonal + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score):
                curr_align1 = [seq1[i - 1]] + align1[::]
                curr_align2 = [seq2[j - 1]] + align2[::]
                paths.append((i - 1, j - 1, curr_align1, curr_align2, curr_path[::]))
            if score_current == score_up + gap_score:
                curr_align1 = [seq1[i - 1]] + align1[::]
                curr_align2 = ["-"] + align2[::]
                paths.append((i - 1, j, curr_align1, curr_align2, curr_path[::]))
            if score_current == score_left + gap_score:
                curr_align1 = ["-"] + align1[::]
                curr_align2 = [seq2[j - 1]] + align2[::]
                paths.append((i, j - 1, curr_align1, curr_align2, curr_path[::]))
        else:
            align_with_paths.append((curr_path, "".join(align1), "".join(align2)))
    return [list(row) for row in score_matrix], max_score, align_with_paths[:2]


def save_result(dp, score, format):
    for row in dp:
        print(row)
    print(score)
    for path, align1, align2 in format:
        print(align1)
        print(align2)
        print(path)


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
    dp, score, paths = find_global_alignments_with_scores(seq1, seq2)
    save_result(dp, score, paths)

    print("------------")
    # Calculate the best end-space free alignment
    dp, score, paths = end_space_free_alignment(seq1, seq2)
    save_result(dp, score, paths)
