import os

current_dir = os.path.dirname(os.path.abspath(__file__))

# Cargar archivos de genomas
def load_genome_file(filename='SARS-COV-2-MT106054.1.txt'):
    file_path = os.path.join(current_dir,'..','data',filename)
    normalized_path = os.path.normpath(file_path)

    genome = {}

    with open(normalized_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                genome_name = line[1:].strip().split('S',1)
                genome_name = genome_name[0]
                genome[genome_name] = ''
            elif not line.startswith('>'):
                genome[genome_name] += line.strip()

    return genome

# Cargar archivos de genes
def load_gene_file(filename):
    file_path = os.path.join(current_dir,'..','data',filename)
    normalized_path = os.path.normpath(file_path)

    gene = {}
    with open(normalized_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                gene_name = line[1:].strip().split('S',1)
                gene_name = gene_name[0]
                gene[gene_name] = ''
            elif not line.startswith('>'):
                gene[gene_name] += line.strip()
    return gene

# Cargar archivo de proteínas
def load_protein_file(filename):
    file_path = os.path.join(current_dir,'..','data',filename)
    normalized_path = os.path.normpath(file_path)

    proteins = {}
    protein_name = None
    with open(normalized_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                protein_name = line[1:].strip()
                proteins[protein_name] = ''
            elif protein_name is not None and not line.startswith('>'):
                proteins[protein_name] += line.strip()

    return proteins

# Limpiar secuencias (espacios, saltos)
def clean_sequence(sequence):
    return 

# Guardar resultados 
def save_results_to_file(data, filename):
    return 