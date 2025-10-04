# PUNTO 2: Palíndromos más largos

'''
Encontrar el palíndromo mas largo en cada uno de los tres genes (gen M, S y ORF1AB). 
    - Los palíndromos son importantes, porque son regiones propensas a mutaciones en un gen. 
Por cada gen, muestra la longitud del palíndromo mas largo y guárdalo en un archivo.
'''

import os
from algorithms.palindrome_finder import manacherAlgorithm

def save_palindrome_to_file(gene_name, palindrome_sequence): 
    results_dir = "results"
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    filename = f"{results_dir}/palindrome_gene_{gene_name}.txt"

    try:
        with open(filename, 'w', encoding="utf-8") as file:
            file.write(
                f"   📏 Palíndromo más largo: {palindrome_sequence['length']}\n"
                f"   📍 Posición en el gen: {palindrome_sequence['position']}\n"
                f"   🔤 Secuencia: {palindrome_sequence['sequence']}\n"
            )

        print(f"   💾 Palíndromo guardado en: {filename}")
        return filename
    
    except Exception as e:
        print(f"   ❌ Error al guardar archivo {filename}: {e}")
        return None

def display_palindrome_results(gene_name, result):
    if result['palindrome_length'] > 0:
        print("-"*50)
        print(f"✅ Gen {gene_name}:")
        print(f"   📏 Palíndromo más largo: {result['palindrome_length']} nucleótidos")
        print(f"   📍 Posición en el gen: {result['palindrome_position']}")
        print(f"   🔤 Secuencia: {result['longest_palindrome'][:50]}{'...' if len(result['longest_palindrome']) > 50 else ''}")
        print(f"   💾 Guardado en: {result['saved_to_file']}")
    else:
        print(f"❌ Gen {gene_name}: No se encontraron palíndromos significativos")

def analyze_palindromes_in_genes(gene_sequence): 
    results = {}
    for gene_name, gene_sequence in gene_sequence.items():
        print(f"\n🧬 Analizando palíndromos en Gen {gene_name}...")

        longest_palindrome = manacherAlgorithm(gene_sequence)
        file_name = save_palindrome_to_file(gene_name,longest_palindrome)

        results[gene_name] = {
            'gene_sequence': gene_sequence,
            'longest_palindrome': longest_palindrome['sequence'],
            'palindrome_length': longest_palindrome['length'],
            'palindrome_position': longest_palindrome['position'],
            'saved_to_file' : file_name
        }

        display_palindrome_results(gene_name,results[gene_name])

    return results

def palindrome_test():
    print("🧪 TEST: Análisis de palíndromos")
    
    # Genes de prueba con palíndromos conocidos
    test_genes = {
        'Test_1': 'ATGCCCGAATTCGGGCATGCATGCCCGAATTCGGGCAT',  # GCCCG
        'Test_2': 'ATGTTTAAATTTCCCGGGGAAATTTAAATTTGGG',      # Contiene TTTAAATTT
        'Test_3': 'ATGGAGCTCGAGCTCGAGCTCGAGCTCCCAT'          # Contiene GAGCTCGAGCTCGAGCTCGAG
    }
    
    # Analizar palíndromos de prueba
    test_results = analyze_palindromes_in_genes(test_genes)