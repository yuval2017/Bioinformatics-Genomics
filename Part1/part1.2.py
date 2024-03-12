import random

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


def generate_random_RNA_sequence(length):
    if length % 5 != 0:
        raise ValueError("Length must be a multiple of 5.")

    # Generate a pool of nucleotides with each nucleotide appearing exactly 5 times
    nucleotides = ['A'] * 5 + ['U'] * 5 + ['C'] * 5 + ['G'] * 5
    random.shuffle(nucleotides)

    # Initialize the RNA sequence
    RNA_sequence = []

    # Generate the RNA sequence ensuring no more than 3 letters in a row are the same
    while len(RNA_sequence) < length:
        random.shuffle(nucleotides)
        RNA_sequence.extend(nucleotides)

    return ''.join(RNA_sequence[:length])


import numpy as np

def nussinov(S, pairs):
    def fill(S, pairs):
        M = len(S)
        dp = np.zeros([M, M], dtype=int)

        for k in range(1, len(S)):
            for i in range(len(S) - k):
                j = i + k
                if j >= i:
                    down = dp[i + 1][j]
                    left = dp[i][j - 1]
                    diag = dp[i + 1][j - 1] + 1 if (S[i], S[j]) in pairs else 0
                    rc = max([dp[i][t] + dp[t + 1][j] for t in range(i, j)])
                    dp[i][j] = max(down, left, diag, rc)
                else:
                    dp[i][j] = 0
        return dp


    def traceback(dp, S, i, L, pairs):
        fold = []
        queue = [(i, L)]
        while queue:
            i, j = queue.pop()
            if i < j:
                if dp[i][j] == dp[i + 1][j]:
                    queue.append((i+1, j))

                elif dp[i][j] == dp[i][j - 1]:
                    queue.append((i, j-1))
                elif dp[i][j] == dp[i + 1][j - 1] + 1 if (S[i], S[j]) in pairs else 0:
                    fold.append((i, j))
                    queue.append((i+1, j-1))

                else:
                    for t in range(i + 1, j - 1):
                        if dp[i][j] == dp[i, t] + dp[t + 1][j]:
                            queue.append((i, t))
                            queue.append((t+1, j))
                            break
        return fold
    def dot_write(S, fold):
        dot = ["."] * len(S)
        for s in fold:
            dot[min(s)] = "("
            dot[max(s)] = ")"
        return "".join(dot)

    dp = fill(S, pairs)
    fold = traceback(dp, S, 0, len(S) - 1, pairs)
    return dp, fold, dot_write(S, fold)




if __name__ == "__main__":
    length = 20
    RNA_sequence = generate_random_RNA_sequence(length)
    print("Random RNA sequence S of length 20:")
    print("The RNA sequence:", RNA_sequence)

    pairs = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}
    table, fold, par = nussinov(RNA_sequence, pairs)

    # Plot the matrix with the optimal secondary structure
    plot_matrix(matrix=table, yellow_squares=fold, square_width=0.055, square_high=0.055)

    print("Optimal RNA secondary structure (base pairs):", fold)
    print("Nussinov algorithm parameters:", par)