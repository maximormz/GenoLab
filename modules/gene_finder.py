# PUNTO 1: Búsqueda de genes en genoma

'''
Encontrar los índices de aparición de cada uno de los tres genes (gen M, S y ORF1AB) en la secuencia del virus SARS-COV-2-MN908947.3.txt. 
Por cada gen muestra su nombre, índices de aparición en la secuencia del virus y primeros 12 caracteres.
'''

from algorithms.string_matching import kmp_search_all_occurrences

def find_genes_index(genome_sequence, gene_sequence): 
    if not genome_sequence:
        print("⚠️  Advertencia: La secuencia del genoma está vacía")
        return {}
    if not gene_sequence:
        print("⚠️  Advertencia: No hay secuencias de genes para buscar")
        return {}
    
    results = {}

    for gene_name,sequence in gene_sequence.items():
        if not sequence:
            print(f"⚠️  Advertencia: Secuencia vacía para el gen {gene_name}")
            continue
        
        direct_indexes = list(kmp_search_all_occurrences(genome_sequence, sequence))

        results[gene_name] = {
            'indexes': direct_indexes,
            'sequence': sequence,
            'first_12': sequence[:12],
            'last_12': sequence[-12:],
            'length': len(sequence)
        }

        display_gene_results(gene_name, results[gene_name])
    
    return results

def display_gene_results(gene_name, result):
    """Muestra los resultados de búsqueda para un gen"""

    if result['indexes']:
        indexes_str = ', '.join(map(str, result['indexes']))
        print(f"✅ Gen {gene_name}: Índices [{indexes_str}], Primeros 12: \"{result['first_12']}...\", Ultimos 12: \",...{result['last_12']}\"")
        print(f"   📏 Longitud del gen: {result['length']} nucleótidos")
        
    else:
        print('-'*50)
        print(f"❌ Gen {gene_name}: No encontrado")
        print(f"   🔍 Buscando: \"{result['first_12']}...{result['last_12']}\"")
        print(f"   📏 Longitud del gen: {result['length']} nucleótidos")

def gene_test():
    print("🧪 TEST: Búsqueda de genes en genoma")
    
    test_genome = "ATGCTAGCTAGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    test_genes = {
        'gen_M': "TAGCTAGATCGA",
        'gen_S': "ATCGATCGATCG", 
        'gen_ORF1AB': "GATCGATCGATC"
    }
    
    print(f"Genoma de prueba: {test_genome[:20]}... (longitud: {len(test_genome)})")
    
    # Ejecutar búsqueda
    results = find_genes_index(test_genome, test_genes)
    
    for gene_name, result in results.items():
        status = "✅ ENCONTRADO" if result['indexes'] else "❌ NO ENCONTRADO"
        print(f"\n{gene_name}: {status}")
        if result['indexes']:
            print(f"   Índices: {result['indexes']}")
            print(f"   Primeros 12: \"{result['first_12']}\"")
            print(f"   Últimos 12: \"{result['last_12']}\"")
        print()