import os

current_dir = os.path.dirname(os.path.abspath(__file__))

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
            if protein_name is not None and not line.startswith('>'):
                proteins[protein_name] += line.strip()

    return proteins