# PUNTO 1: Búsqueda de genes en genoma

'''
Encontrar los índices de aparición de cada uno de los tres genes (gen M, S y ORF1AB) en la secuencia del virus SARS-COV-2-MN908947.3.txt. 
Por cada gen muestra su nombre, índices de aparición en la secuencia del virus y primeros 12 caracteres.
'''

from algorithms.string_matching import kmp_search_all_occurrences

def find_genes_index(genome_sequence, gene_sequence): 
    for gene_name in gene_sequence.keys():
        Index = list(kmp_search_all_occurrences(genome_sequence['MT106054.1'], gene_sequence[gene_name]))
        print(f"Gen {gene_name}: Índices {Index}, Primeros 12: \"{gene_sequence[gene_name][:12]}\"")