# DecodingDNA

**DecodingDNA**
- This project offers a set of utility functions designed to perform common tasks in bioinformatics such as sequence alignment, finding coding sequences, analyzing GC content, and codon usage. It utilizes the Biopython library and matplotlib for visualization.

**Key Features**
- Sequence Alignment: The sequence_alignment function performs pairwise sequence alignment using the Biopython's pairwise2 module. It uses the globalxx alignment method which finds the alignment with the highest similarity (maximum number of matches and minimum number of mismatches).
- Coding Sequence Extraction: The find_coding_sequence function identifies the coding sequence in a given DNA sequence. It locates the start codon (ATG) and any of the three stop codons (TAA, TAG, TGA). It returns the coding sequence or an error message if either the start or stop codon can't be found.
- GC Content Analysis: The measure_gc function calculates the GC content of a DNA sequence, optionally considering only the coding sequence. It can be used to get a quick overview of the composition of a DNA sequence.

- Codon Usage Analysis: The analyze_codon_usage function provides an analysis of codon usage in a DNA sequence. It returns a dictionary where the keys are codons and the values are their relative frequencies.

**Visualization**
- Two visualization functions are included:
  - plot_gc_content_distribution: This function takes a list of DNA sequences and plots a histogram of GC content distribution. It can consider either the entire sequence or only the coding sequence based on the coding parameter.
  - plot_codon_usage_bias: This function takes a list of DNA sequences and plots a heatmap showing the relative frequencies of each codon across all sequences. It uses matplotlib to create an interactive heatmap with color coding to show frequency distribution.

**Usage**
- Each function in the module can be used independently. Simply import the desired functions into your script and provide the necessary arguments. Most functions require a DNA sequence in string format.
- This toolkit provides fundamental utilities that can be expanded or incorporated into more complex bioinformatics analyses. It serves as a great starting point for those new to bioinformatics or as a quick reference for more seasoned researchers.
