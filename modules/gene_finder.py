# PUNTO 1: Búsqueda de genes en genoma

'''
Encontrar los índices de aparición de cada uno de los tres genes (gen M, S y ORF1AB) en la secuencia del virus SARS-COV-2-MN908947.3.txt. 
Por cada gen muestra su nombre, índices de aparición en la secuencia del virus y primeros 12 caracteres.
'''

from algorithms.string_matching import kmp_search_all_occurrences

def get_reverse_complement(sequence):
    complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    
    # Crear complemento
    complement = ''.join(complement_map.get(base, base) for base in sequence.upper())
    
    # Reverso del complemento
    reverse_complement = complement[::-1]
    
    return reverse_complement

def find_genes_index(genome_sequence, gene_sequence): 
    genome = genome_sequence['MN908947.3']
    results = {}

    for gene_name in gene_sequence.keys():
        
        Index = list(kmp_search_all_occurrences(genome, gene_sequence[gene_name]))

        reverse_comp_index = []
        if not Index:
            reverse_complement = get_reverse_complement(gene_sequence[gene_name])
            reverse_comp_index = list(kmp_search_all_occurrences(genome, reverse_complement))
        
        all_indexes = Index + reverse_comp_index
        results[gene_name] = {
            'indexes': all_indexes,
            'sequence': gene_sequence[gene_name],
            'first_12': gene_sequence[gene_name][:12],
            'found_as': 'indexes' if Index else 'reverse_complement' if reverse_comp_index else 'not found'
        }

        display_gene_results(gene_name, results[gene_name])
    
    return results

def display_gene_results(gene_name, result):
    """Muestra los resultados de búsqueda para un gen"""
    if result['indexes']:
        indexes_str = ', '.join(map(str, result['indexes']))
        print(f"✅ Gen {gene_name}: Índices [{indexes_str}], Primeros 12: \"{result['first_12']}\"")
        
        if result['found_as'] == 'reverse_complement':
            print(f"   ⚠️  Encontrado como complemento reverso")
    else:
        print(f"❌ Gen {gene_name}: No encontrado")
        print(f"   🔍 Buscando: \"{result['first_12']}...\"")
        print(f"   📏 Longitud del gen: {len(result['sequence'])} nucleótidos")