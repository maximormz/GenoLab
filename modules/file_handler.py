import os

current_dir = os.path.dirname(os.path.abspath(__file__))

# Cargar archivos de genomas
def load_genome_file(filename):
    file_path = os.path.join(current_dir,'..','data',filename)
    normalized_path = os.path.normpath(file_path)

    genome = {}

    try:
        with open(normalized_path, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    genome_name = line[1:].strip().split('S',1)
                    genome_name = genome_name[0]
                    genome[genome_name] = ''
                elif not line.startswith('>'):
                    genome[genome_name] += line.strip()

        return genome
    except FileNotFoundError:
        print(f"Error: El archivo '{filename}' no se encontró en {normalized_path}.")
        return {}

# Cargar archivos de genes
def load_gene_file(filename):
    file_path = os.path.join(current_dir,'..','data',filename)
    normalized_path = os.path.normpath(file_path)

    gene = {}

    try:
        with open(normalized_path, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    gene_name = line[1:].strip().split('S',1)
                    gene_name = gene_name[0]
                    gene[gene_name] = ''
                elif not line.startswith('>'):
                    gene[gene_name] += line.strip()
        return gene
    except FileNotFoundError:
        print(f"Error: El archivo '{filename}' no se encontró en {normalized_path}.")
        return {}

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
        return {}

# Limpiar secuencias (espacios, saltos)
def clean_sequence(sequence):
    return 

# Guardar resultados 
def save_results_to_file(data, filename):
    return 