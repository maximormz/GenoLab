from modules.file_handler import load_all_files
from modules.gene_finder import find_genes_index

if __name__ == "__main__":
    print("===Analisis Genomico de SARS-COV-2===\n")

    print("Cargando archivos...")
    files_data = load_all_files()
    if files_data:
        print("Archivos cargados exitosamente!\n")

    print("---PUNTO 1: Busqueda de genes en genoma---")
    
    find_genes_index(files_data['genome'], files_data['genes'])