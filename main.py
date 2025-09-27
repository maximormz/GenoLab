from modules.file_handler import load_all_files
from modules.gene_finder import gene_test, find_genes_index
from modules.palindrome_analyzer import palindrome_test, analyze_palindromes_in_genes

if __name__ == "__main__":
    print("===Analisis Genomico de SARS-COV-2===\n")

    print("Cargando archivos...")
    files_data = load_all_files()
    if files_data:
        print("Archivos cargados exitosamente!\n")
        print("-"*50)

    print("🎯 PUNTO 1: Búsqueda de Genes")
    print("=" * 50)

    #gene_test()
    #find_genes_index(files_data['genome']['MN908947.3'], files_data['genes'])

    print("\n🔄 PUNTO 2: Palíndromos")
    print("=" * 50)

    #palindrome_test()
    #analyze_palindromes_in_genes(files_data['genes'])

    print("\n🧬 PUNTO 3: Mapeo de Proteínas")
    print("=" * 50)

    print("\n⚖️ PUNTO 4: Comparación de Genomas")
    print("=" * 50)