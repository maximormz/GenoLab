import os

current_dir = os.path.dirname(os.path.abspath(__file__))

# Cargar archivos de genomas
def load_genome_file(filename):
    return 

# Cargar archivos de genes
def load_gene_file(filename='gen-M.txt'):
    file_path = os.path.join(current_dir,'..','data',filename)
    normalized_path = os.path.normpath(file_path)

    gene = {}
    gene_name = None
    with open(normalized_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                gene_name = line[1:].strip().split('S',1)
                gene_name = gene_name[0]
                gene[gene_name] = ''
            elif gene_name is not None and not line.startswith('>'):
                gene[gene_name] += line.strip()
    return gene

GEN = load_gene_file()
for key in GEN:
    print(F'{key}: {GEN[key]}')

'''
Encontrar el palíndromo mas largo en cada uno de los tres genes (gen M, S y ORF1AB). Los palíndromos son importantes, porque son regiones propensas a mutaciones en un gen. 

Por cada gen, muestra la longitud del palíndromo mas largo y guárdalo en un archivo.
'''

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