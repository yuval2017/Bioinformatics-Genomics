import random
import numpy as np
from matplotlib import pyplot as plt


def plot_matrix(matrix, square_width=0.1, square_high=0.1, yellow_squares=[], green_squares=[]):
    fig, ax = plt.subplots()
    ax.axis('tight')
    ax.axis('off')
    the_table = ax.table(cellText=matrix, loc='center')

    # Increase cell padding to make each square bigger
    cell_dict = the_table.get_celld()
    for key in cell_dict:
        cell = cell_dict[key]
        cell.set_height(square_high)  # Adjust this value as needed
        cell.set_width(square_width)  # Adjust this value as needed
        cell.get_text().set_fontsize(16)  # Adjust this value as needed for the font size

        # Center the text horizontally and vertically within each cell
        cell.get_text().set_horizontalalignment('center')
        cell.get_text().set_verticalalignment('center')

    # Paint yellow and green squares
    for square in yellow_squares:
        if square in green_squares:
            cell = cell_dict[square]
            cell.set_facecolor('yellowgreen')
        else:
            cell = cell_dict[square]
            cell.set_facecolor('yellow')

    for square in green_squares:
        if square not in yellow_squares:
            cell = cell_dict[square]
            cell.set_facecolor('green')

    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()


def process_and_plot(func, seq1, seq2):
    dp, score, paths = func(seq1, seq2)
    print(paths[0][1])
    path1, str_seq11, str_seq12 = (paths[0][0], paths[0][1], paths[0][2]) if paths else ([], "", "")
    path2, str_seq21, str_seq22 = (paths[1][0], paths[1][1], paths[1][2]) if len(paths) > 1 else ([], "", "")
    plot_matrix(matrix=dp, yellow_squares=path1, green_squares=path2)
    print("Traces:")
    print("Trace1:", path1)
    print("Trace2:", path2)
    print("First trace alignments:")
    print("Alignment1:", str_seq11)
    print("Alignment1:", str_seq12)
    print("Second trace alignments:")
    print("Alignment1:", str_seq21)
    print("Alignment1:", str_seq22)

    print("Score:", score)


# Function to generate a random DNA sequence of given length
import random


def generate_dna_sequences(length=10):
    # Define the DNA alphabet
    dna_alphabet = ['A', 'C', 'T', 'G']

    # Generate the first DNA sequence ensuring each nucleotide appears at least once
    s1 = random.sample(dna_alphabet, 4)  # Take one of each nucleotide
    remaining_letters = random.choices(dna_alphabet, k=length - 4)  # Fill the rest randomly
    random.shuffle(remaining_letters)  # Shuffle the remaining letters
    s1 += remaining_letters

    # Generate the second DNA sequence ensuring each nucleotide appears at least once
    s2 = random.sample(dna_alphabet, 4)  # Take one of each nucleotide
    remaining_letters = random.choices(dna_alphabet, k=length - 4)  # Fill the rest randomly
    random.shuffle(remaining_letters)  # Shuffle the remaining letters
    s2 += remaining_letters

    # Shuffle both sequences to randomize the order
    random.shuffle(s1)
    random.shuffle(s2)

    # Convert lists to strings
    s1 = ''.join(s1)
    s2 = ''.join(s2)

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
        trace.append((0, 0))
        format.append((trace, align1, align2))

    return dp, dp[n][m], format


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


# Example usage

# Example usage
# Test the function with example sequences


if __name__ == "__main__":
    # Generate sequences of length 10
    seq1, seq2 = generate_dna_sequences()
    print("First DNA sequence (S1):", seq1)
    print("Second DNA sequence (S2):", seq2)
    process_and_plot(find_global_alignments_with_scores, seq1, seq2)
    print("------------")
    process_and_plot(end_space_free_alignment, seq1, seq2)
