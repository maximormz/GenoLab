import os

current_dir = os.path.dirname(os.path.abspath(__file__))

# Cargar archivos de genomas
def load_genome_file(filename):
    file_path = os.path.join(current_dir,'..','data',filename)
    normalized_path = os.path.normpath(file_path)

    genome = ''

    try:
        with open(normalized_path, 'r') as file:
            for line in file:
                if not line.startswith('>'):
                    genome += line.strip("\n")

        return genome
    except FileNotFoundError:
        print(f"Error: El archivo '{filename}' no se encontró en {normalized_path}.")
        return ""

# Cargar archivos de genes
def load_gene_file(filename):
    file_path = os.path.join(current_dir,'..','data',filename)
    normalized_path = os.path.normpath(file_path)
    gene = ''
    try:
        with open(normalized_path, 'r') as file:
            for line in file:
                if not line.startswith('>'):
                    gene += line.strip("\n")
        return gene
    except FileNotFoundError:
        print(f"Error: El archivo '{filename}' no se encontró en {normalized_path}.")
        return ""

# Cargar archivo de proteínas
def load_protein_file(filename):
    file_path = os.path.join(current_dir,'..','data',filename)
    normalized_path = os.path.normpath(file_path)

    proteins = {}
    protein_name = None

    try:
        with open(normalized_path, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    protein_name = line[1:].strip()
                    proteins[protein_name] = ''
                elif protein_name is not None and not line.startswith('>'):
                    proteins[protein_name] += line.strip()

        return proteins
    except FileNotFoundError:
        print(f"Error: El archivo '{filename}' no se encontró en {normalized_path}.")
        return ""

# Retornar todos los archivos cargados
def load_all_files():
    genome_data_1 = load_genome_file('SARS-COV-2-MN908947.3.txt')
    genome_data_2 = load_genome_file('SARS-COV-2-MT106054.1.txt')
    gene_M_data = load_gene_file('gen-M.txt')
    gene_S_data = load_gene_file('gen-S.txt')
    gene_ORF1AB_data = load_gene_file('gen-ORF1AB.txt')
    protein_data = load_protein_file('seq-proteins.txt')

    return {
        'genome': {
            'MN908947.3' : genome_data_1,
            'MT106054.1' : genome_data_2},
        'genes': {
            'M': gene_M_data,
            'S': gene_S_data,
            'ORF1AB': gene_ORF1AB_data
        },
        'proteins': protein_data
    }

# Guardar resultados 
def save_results_to_file(data, filename):
    return 