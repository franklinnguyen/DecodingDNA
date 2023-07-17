import numpy as np
import matplotlib.pyplot as plt
from Bio import pairwise2


def sequence_alignment(sequence1, sequence2):
    """
    Performs a pairwise sequence alignment between two sequences.

    Args:
        sequence1 (str): The first sequence to align.
        sequence2 (str): The second sequence to align.

    Returns:
        None: This function does not return a value.

    Raises:
        None.

    This function uses the Biopython library to perform a pairwise sequence alignment
    between two input sequences. It applies the global alignment algorithm using the
    'x' (matches/mismatches) substitution matrix.

    The function prints the alignments between the two sequences, including the aligned
    sequences and the alignment score.

    Example:
        sequence1 = "ACGTCGTACG"
        sequence2 = "ACGGTACGTT"
        sequence_alignment(sequence1, sequence2)
    """
    alignments = pairwise2.align.globalxx(sequence1, sequence2)
    # Print the alignments
    for alignment in alignments:
        print(pairwise2.format_alignment(*alignment))


def find_coding_sequence(original):
    """
    Finds the coding sequence in the given DNA sequence.

    Args:
        original (str): The original DNA sequence.

    Returns:
        str: The coding sequence if found.
        str: 'ERROR, start codon not found!' if the start codon is not found.
        str: 'ERROR, stop codon not found!' if the stop codon is not found.
    """
    start = "ATG"
    stop = {"TAA", "TAG", "TGA"}

    start_index = original.find(start)

    if start_index == -1:
        return "ERROR, start codon not found!"

    # sequence beginning at start
    shortened_sequence = original[start_index:]
    triplets = split_triplets(shortened_sequence)

    stop_index = None

    for i, codon in enumerate(triplets):
        if codon in stop:
            stop_index = i
            break

    if stop_index is None:
        return "ERROR, stop codon not found!"

    coding_sequence = "".join(shortened_sequence[:stop_index])

    return coding_sequence


def split_triplets(sequence):
    """
    Splits the given sequence into triplets (codons).

    Args:
        sequence (str): The input sequence.

    Returns:
        list: A list of triplets (codons).
    """
    return [sequence[i : i + 3] for i in range(0, len(sequence), 3)]


def measure_gc(sequence, coding=True):
    """
    Calculates the GC content (percentage of G and C bases) in the given DNA sequence.

    Args:
        sequence (str): The DNA sequence to analyze.
        coding (bool, optional): Indicates whether to consider only the coding sequence.
                                 Defaults to True.

    Returns:
        float: The GC content as a ratio.
    """
    if coding:
        sequence = find_coding_sequence(sequence)
    total_bases = len(sequence)
    gc_content = sequence.count("G") + sequence.count("C")
    return gc_content / total_bases


def analyze_codon_usage(sequence):
    """
    Analyzes the codon usage in the given DNA sequence.

    Args:
        sequence (str): The DNA sequence to analyze.

    Returns:
        dict: A dictionary where the keys are codons and the values are their frequencies.
    """
    triplets = split_triplets(sequence)
    codon_count = {}
    total_codons = len(triplets)

    for triplet in triplets:
        if triplet not in codon_count:
            codon_count[triplet] = 1
        else:
            codon_count[triplet] += 1

    codon_frequencies = {key: val / total_codons for key, val in codon_count.items()}

    return codon_frequencies


def plot_gc_content_distribution(sequences, coding):
    """
    Plots the distribution of GC content among the given DNA sequences.

    Args:
        sequences (list): A list of DNA sequences.
        coding (bool, optional): Indicates whether to consider only the coding sequence
                                 when calculating GC content. Defaults to True.

    Returns:
        None: This function displays the histogram plot.

    Raises:
        None.

    The function calculates the GC content for each sequence and plots a histogram
    showing the distribution of GC content percentages.

    The x-axis represents the GC content (%) ranging from 0 to 100. The y-axis represents
    the frequency of sequences falling within each GC content bin.

    The `sequences` parameter should be a list of DNA sequences in string format.

    If `coding` is set to True (default), the function calculates the GC content
    considering only the coding sequence (between the start codon and the stop codon).
    If `coding` is False, the GC content is calculated for the entire sequence.

    Example:
        sequences = ['ATCGATCG', 'GCGCTAGC', 'ATATATAT']
        plot_gc_content_distribution(sequences, coding=True)
    """
    gc_contents = [measure_gc(sequence, coding) for sequence in sequences]

    plt.hist(gc_contents, bins=20, edgecolor="black", color="darkseagreen")
    plt.xlabel("GC Content (%)")
    plt.ylabel("Frequency")
    plt.title("GC Content Distribution")
    plt.show()


import matplotlib.pyplot as plt
from collections import OrderedDict


def plot_codon_usage_bias(sequences):
    """
    Plots the codon usage bias among the given DNA sequences.

    Args:
        sequences (list): A list of DNA sequences.

    Returns:
        None: This function displays the bar plot.

    Raises:
        None.

    The function analyzes the codon usage in each sequence and plots a bar graph
    showing the relative frequencies of each codon.

    The x-axis represents the codons, and the y-axis represents the frequency of
    each codon.

    The `sequences` parameter should be a list of DNA sequences in string format.

    Example:
        sequences = ['ATCGATCG', 'GCGCTAGC', 'ATATATAT']
        plot_codon_usage_bias(sequences)
    """
    codon_frequencies = [analyze_codon_usage(sequence) for sequence in sequences]
    codons = list(codon_frequencies[0].keys())

    frequencies = np.array(
        [[freq.get(codon, 0) for codon in codons] for freq in codon_frequencies]
    )
    frequencies = frequencies.reshape((len(sequences), -1))

    fig, ax = plt.subplots()
    im = ax.imshow(frequencies, cmap="Blues")

    codon_frequencies = [analyze_codon_usage(sequence) for sequence in sequences]
    codons = list(codon_frequencies[0].keys())

    frequencies = np.array(
        [[freq.get(codon, 0) for codon in codons] for freq in codon_frequencies]
    )
    frequencies = frequencies.reshape((len(sequences), -1))

    fig, ax = plt.subplots()
    im = ax.imshow(frequencies, cmap="Blues")

    # Set ticks and labels for x-axis
    ax.set_xticks(np.arange(len(codons)))
    ax.set_xticklabels(codons)

    # Set ticks and labels for y-axis
    ax.set_yticks(np.arange(len(sequences)))
    ax.set_yticklabels(np.arange(1, len(sequences) + 1))

    # Rotate the x-axis labels for better visibility
    plt.xticks(rotation=90)

    # Add colorbar
    cbar = ax.figure.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
    cbar.ax.set_ylabel("Frequency", rotation=-90, va="bottom")

    # Set plot title and axis labels
    plt.title("Codon Usage Bias")
    plt.xlabel("Codons")
    plt.ylabel("Sequences")

    # Display the plot
    plt.show()


if __name__ == "__main__":
    pass