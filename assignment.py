from Bio.Seq import Seq
from collections import Counter


def dna_to_amino_acids(dna_sequence):
    
    dna_seq = Seq(dna_sequence)
    
    mRNA_seq = dna_seq.transcribe()
    
    amino_acids = mRNA_seq.translate()
    
    return mRNA_seq, amino_acids

dna_sequence = "TTACGA"
mRNA_seq, amino_acids = dna_to_amino_acids(dna_sequence)

print(f"Input DNA: {dna_sequence}")
print(f"mRNA Sequence: {mRNA_seq}")
print(f"Amino Acids: {amino_acids}")

codon_table = {
    'F': ['UUU', 'UUC'], 'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    'I': ['AUU', 'AUC', 'AUA'], 'M': ['AUG'], 'V': ['GUU', 'GUC', 'GUA', 'GUG'],
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'], 'P': ['CCU', 'CCC', 'CCA', 'CCG'],
    'T': ['ACU', 'ACC', 'ACA', 'ACG'], 'A': ['GCU', 'GCC', 'GCA', 'GCG'],
    'Y': ['UAU', 'UAC'], 'H': ['CAU', 'CAC'], 'Q': ['CAA', 'CAG'], 'N': ['AAU', 'AAC'],
    'K': ['AAA', 'AAG'], 'D': ['GAU', 'GAC'], 'E': ['GAA', 'GAG'], 'C': ['UGU', 'UGC'],
    'W': ['UGG'], 'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'G': ['GGU', 'GGC', 'GGA', 'GGG'],
    'STOP': ['UAA', 'UAG', 'UGA']
}

def get_codon_frequencies(amino_acid_sequence, mRNA_sequence):
    
    possible_codons = []
    for amino_acid in amino_acid_sequence:
        codons = codon_table[amino_acid]
        possible_codons.extend(codons)
    
    codons_in_mRNA = [mRNA_sequence[i:i+3] for i in range(0, len(mRNA_sequence), 3)]
    
    codon_count = Counter(codons_in_mRNA)
    
    relevant_codon_count = {codon: codon_count[codon] for codon in possible_codons if codon in codon_count}
    
    return relevant_codon_count

amino_acid_sequence = "WYW"
mRNA_sequence = "UGGGUACUGG"

codon_frequencies = get_codon_frequencies(amino_acid_sequence, mRNA_sequence)

print(f"\nAmino Acid Sequence: {amino_acid_sequence}")
print(f"mRNA Sequence: {mRNA_sequence}")
print("Codon Frequencies:")
for codon, freq in codon_frequencies.items():
    print(f"{codon}: {freq}")
