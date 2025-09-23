
GENETIC_CODE = {
    # Fenilalanina (F)
    'TTT': 'F', 'TTC': 'F',
    
    # Leucina (L)
    'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    
    # Serina (S)
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
    
    # Tirosina (Y)
    'TAT': 'Y', 'TAC': 'Y',
    
    # Codones de parada (Stop)
    'TAA': '*', 'TAG': '*', 'TGA': '*',
    
    # Cisteína (C)
    'TGT': 'C', 'TGC': 'C',
    
    # Triptófano (W)
    'TGG': 'W',
    
    # Prolina (P)
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    
    # Histidina (H)
    'CAT': 'H', 'CAC': 'H',
    
    # Glutamina (Q)
    'CAA': 'Q', 'CAG': 'Q',
    
    # Arginina (R)
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    
    # Isoleucina (I)
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
    
    # Metionina (M) - Codón de inicio
    'ATG': 'M',
    
    # Treonina (T)
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    
    # Asparagina (N)
    'AAT': 'N', 'AAC': 'N',
    
    # Lisina (K)
    'AAA': 'K', 'AAG': 'K',
    
    # Valina (V)
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    
    # Alanina (A)
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    
    # Ácido aspártico (D)
    'GAT': 'D', 'GAC': 'D',
    
    # Ácido glutámico (E)
    'GAA': 'E', 'GAG': 'E',
    
    # Glicina (G)
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}                                        # Tabla completa de codones

START_CODONS = ['ATG']                   # Codones de inicio

STOP_CODONS = ['TAA', 'TAG', 'TGA']      # Codones de parada

AMINO_ACID_TO_CODONS = {}
for codon, amino_acid in GENETIC_CODE.items():
    if amino_acid not in AMINO_ACID_TO_CODONS:
        AMINO_ACID_TO_CODONS[amino_acid] = []
    AMINO_ACID_TO_CODONS[amino_acid].append(codon)

AMINO_ACIDS = {
    'A': 'Alanina',
    'R': 'Arginina', 
    'N': 'Asparagina',
    'D': 'Ácido aspártico',
    'C': 'Cisteína',
    'Q': 'Glutamina',
    'E': 'Ácido glutámico',
    'G': 'Glicina',
    'H': 'Histidina',
    'I': 'Isoleucina',
    'L': 'Leucina',
    'K': 'Lisina',
    'M': 'Metionina',
    'F': 'Fenilalanina',
    'P': 'Prolina',
    'S': 'Serina',
    'T': 'Treonina',
    'W': 'Triptófano',
    'Y': 'Tirosina',
    'V': 'Valina',
    '*': 'Codón de parada'
}

def codon_to_amino_acid(codon):
    """Convierte un codón de ADN a su aminoácido correspondiente."""
    codon = codon.upper()
    return GENETIC_CODE.get(codon, f'No se ha encontrado ningun aminoacido con esa referencia: {codon}')

def is_start_codon(codon):
    """Verifica si un codón es de inicio"""
    return codon.upper() in START_CODONS

def is_stop_codon(codon):
    """Verifica si un codón es de parada"""
    return codon.upper() in STOP_CODONS