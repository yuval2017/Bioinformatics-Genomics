import random


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
    print(len(RNA_sequence), RNA_sequence)
    # Example usage:
    sequence = "GGGAAAUCC"
    sequence2 = "GGGAAAUCC"

    pairs = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}
    ell = 1
    table, fold, par = nussinov(sequence, pairs)
    for row in table:
        print(row)

    print("Optimal RNA secondary structure (base pairs):", fold)
    print(par)