from collections import Counter

def dna_to_mrna(dna_sequence):
    return dna_sequence.replace('T', 'U')

codon_to_amino_acid = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'UAU': 'Y', 'UAC': 'Y', 'CAU': 'H', 'CAC': 'H',
    'CAA': 'Q', 'CAG': 'Q', 'AAU': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'GAU': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E', 'UGU': 'C', 'UGC': 'C',
    'UGG': 'W', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R',
    'CGG': 'R', 'AGU': 'S', 'AGC': 'S', 'AGA': 'R',
    'AGG': 'R', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G', 'UAA': 'STOP', 'UAG': 'STOP', 'UGA': 'STOP'
}

def translate_mrna_to_amino_acids(mrna_sequence):
    amino_acids = []
    for i in range(0, len(mrna_sequence), 3):
        codon = mrna_sequence[i:i+3]
        amino_acid = codon_to_amino_acid.get(codon, '?')  
        if amino_acid == 'STOP':
            break  
        amino_acids.append(amino_acid)
    return ''.join(amino_acids)

def get_codon_frequencies(amino_acid_sequence, mRNA_sequence):
    
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

    codons_in_mRNA = [mRNA_sequence[i:i+3] for i in range(0, len(mRNA_sequence), 3)]

    codon_count = Counter(codons_in_mRNA)

    relevant_codon_count = {}
    for amino_acid in amino_acid_sequence:
        codons = codon_table.get(amino_acid, [])
        for codon in codons:
            if codon in codon_count:
                relevant_codon_count[codon] = codon_count[codon]

    return relevant_codon_count

dna_sequence = "TTACGA"
mRNA_sequence = dna_to_mrna(dna_sequence)
amino_acids = translate_mrna_to_amino_acids(mRNA_sequence)

print(f"DNA: {dna_sequence}")
print(f"mRNA: {mRNA_sequence}")
print(f"Amino Acids: {amino_acids}")

amino_acid_sequence = "WYW"
mRNA_sequence = "UGGGUACUGG"
codon_frequencies = get_codon_frequencies(amino_acid_sequence, mRNA_sequence)

print(f"\nAmino Acid Sequence: {amino_acid_sequence}")
print(f"mRNA Sequence: {mRNA_sequence}")
print("Codon Frequencies:")
for codon, freq in codon_frequencies.items():
    print(f"{codon}: {freq}")
