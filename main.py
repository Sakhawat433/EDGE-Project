import re

def calculate_gc_content(dna_sequence):
    """
    Calculate the GC content of a DNA sequence.

    :param dna_sequence: str, DNA sequence (e.g., "ATGC")
    :return: float, GC content percentage
    """
    gc_count = sum(1 for base in dna_sequence.upper() if base in "GC")
    return (gc_count / len(dna_sequence)) * 100

def transcribe_dna_to_rna(dna_sequence):
    """
    Transcribe a DNA sequence to RNA.

    :param dna_sequence: str, DNA sequence (e.g., "ATGC")
    :return: str, RNA sequence (e.g., "AUGC")
    """
    return dna_sequence.upper().replace("T", "U")

def find_orfs(dna_sequence, min_length=100):
    """
    Find all open reading frames (ORFs) in a DNA sequence.

    :param dna_sequence: str, DNA sequence (e.g., "ATG...")
    :param min_length: int, minimum ORF length in base pairs
    :return: list of tuples, each tuple contains (start_position, end_position, orf_sequence)
    """
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []

    # Search for ORFs in the 3 reading frames
    for frame in range(3):
        sequence = dna_sequence[frame:]
        for match in re.finditer(f"{start_codon}.*?({'|'.join(stop_codons)})", sequence):
            start = match.start() + frame
            end = match.end() + frame
            orf = dna_sequence[start:end]
            if len(orf) >= min_length:
                orfs.append((start, end, orf))

    return orfs

# Example usage
def main():
    dna_sequence = "ATGCGTATAGCGCTTAAATGCGCTGA"
    print("DNA Sequence:", dna_sequence)

    # Calculate GC Content
    gc_content = calculate_gc_content(dna_sequence)
    print(f"GC Content: {gc_content:.2f}%")

    # Transcribe DNA to RNA
    rna_sequence = transcribe_dna_to_rna(dna_sequence)
    print("RNA Sequence:", rna_sequence)

    # Find ORFs
    orfs = find_orfs(dna_sequence, min_length=9)
    print("Open Reading Frames:")
    for start, end, orf in orfs:
        print(f"Start: {start}, End: {end}, ORF: {orf}")

if __name__ == "__main__":
    main()
