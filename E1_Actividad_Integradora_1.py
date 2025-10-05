import os
import re

# CONSTANTES
current_dir = os.path.dirname(os.path.abspath(__file__))

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}                                        # Tabla completa de codones

START_CODONS = ['ATG']                   # Codones de inicio

STOP_CODONS = ['TAA', 'TAG', 'TGA']      # Codones de parada

# ------------------------------------------------------------------------

# LECTURA AUTOMATICA DE ARCHIVOS
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

# ------------------------------------------------------------------------

# FUNCIONES DE APOYO
def kmp_search_all_occurrences(text, pattern): # Retorna un generador (Convertir en lista al recibirlo)
    if not text or not pattern:
        return
    
    n = len(text)
    m = len(pattern)

    i,j = 0,0

    failure_table = build_failure_table(pattern)

    while i < n:
        if text[i] == pattern[j]:
            i += 1
            j += 1

            if j == m:
                yield i - m
                j = failure_table[j-1]
        else:
            if j > 0:
                j = failure_table[j-1]
            else:
                i += 1

def build_failure_table(pattern):
    """Construye tabla de fallos para KMP"""
    n = len(pattern)
    lps = [0] * n

    length = 0
    i = 1

    while i < n:
        if pattern[i] == pattern[length]:
            length += 1
            lps[i] = length
            i += 1
        else:
            if length != 0:
                length = lps[length - 1]
            else:
                lps[i] = 0
                i += 1
    return lps

def new_phrase(text):
    n = len(text)
    e = 2 * n + 3

    new_string = ["$"] * e
    new_string[0] = "@"
    new_string[-1] = "#"

    for i in range(n):
        new_string[2 * i + 2] = text[i]

    return e, new_string

def manacherAlgorithm(text):
    centro = limite = 0
    e,string = new_phrase(text)
    P = [0] * e

    for i in range(1,e - 1):
        if i < limite:
            simetrica = 2 * centro - i
            P[i] = min(limite - i,P[simetrica])

        gap = P[i] + 1
        while string[i-gap] == string[i+gap]:
            P[i] += 1
            gap += 1
        
        if i + P[i] > limite:
            limite = i + P[i]
            centro = i

    max_len = max(P)
    center_index = P.index(max_len)

    start = (center_index - max_len) // 2
    end = start + max_len

    longest_palindrome = text[start:end]

    return {
        'sequence': longest_palindrome,
        'length': max_len,
        'position': [start,end-1],
    }

# ------------------------------------------------------------------------

# FUNCIONES DEL PUNTO 1
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
        print(f"✅ Gen {gene_name}: Índices [{indexes_str}], \n   🔤 Primeros y ultimos 12: \"{result['first_12']}...{result['last_12']}\"")
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

# ------------------------------------------------------------------------

# FUNCIONES DEL PUNTO 2
def save_palindrome_to_file(gene_name, palindrome_sequence): 
    results_dir = "results"
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    filename = f"{results_dir}/palindrome_gene_{gene_name}.txt"

    try:
        with open(filename, 'w', encoding="utf-8") as file:
            file.write(
                f"   📏 Palíndromo más largo: {palindrome_sequence['length']}\n"
                f"   📍 Posición del Gen | inicial: {palindrome_sequence['position'][0]} final: {palindrome_sequence['position'][1]} \n"
                f"   🔤 Secuencia: {palindrome_sequence['sequence']}\n"
            )
        return filename
    
    except Exception as e:
        print(f"   ❌ Error al guardar archivo {filename}: {e}")
        return None

def display_palindrome_results(gene_name, result):
    if result['palindrome_length'] > 0:
        print("-"*50)
        print(f"✅ Gen {gene_name}:")
        print(f"   📏 Palíndromo más largo: {result['palindrome_length']} nucleótidos")
        print(f"   📍 Posición del Gen | inicial: {result['palindrome_position'][0]} final: {result['palindrome_position'][1]}")
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

# ------------------------------------------------------------------------

# FUNCIONES DEL PUNTO 3
def read_fasta(path):
    header = None
    seq_lines = []
    with open(path, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                header = line[1:].strip()
            else:
                seq_lines.append(line.strip())
    seq = ''.join(seq_lines).upper()
    seq = re.sub('[^ATCGN]', '', seq)
    return header, seq

def read_proteins(path):
    proteins = {}
    cur_id = None
    cur_lines = []
    with open(path, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if cur_id:
                    seq = ''.join(cur_lines).upper()
                    seq = re.sub('[^A-Z]', '', seq)
                    proteins[cur_id] = seq
                cur_id = line.lstrip('>').split()[0]
                cur_lines = []
            else:
                cur_lines.append(line.strip())
    if cur_id:
        seq = ''.join(cur_lines).upper()
        seq = re.sub('[^A-Z]', '', seq)
        proteins[cur_id] = seq
    return proteins

def translate_dna(dna, frame=0):
    aa = []
    for i in range(frame, len(dna)-2, 3):
        codon = dna[i:i+3]
        aa.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(aa)

def kmp_search_all(text, pattern):
    if not pattern:
        return []
    # Construccion de LPS
    m = len(pattern)
    lps = [0]*m
    length = 0
    i = 1
    while i < m:
        if pattern[i] == pattern[length]:
            length += 1
            lps[i] = length
            i += 1
        else:
            if length != 0:
                length = lps[length-1]
            else:
                lps[i] = 0
                i += 1
    # Search
    res = []
    n = len(text)
    i = j = 0
    while i < n:
        if text[i] == pattern[j]:
            i += 1
            j += 1
            if j == m:
                res.append(i-j)
                j = lps[j-1]
        else:
            if j != 0:
                j = lps[j-1]
            else:
                i += 1
    return res

def map_protein_to_genome(protein_seq, genome_seq):
    protein_seq = protein_seq.upper()
    results = []
    for frame in range(3):
        translated = translate_dna(genome_seq, frame)
        occs = kmp_search_all(translated, protein_seq)
        for occ in occs:
            start_nt = frame + occ*3 + 1    # 1-based
            end_nt = start_nt + len(protein_seq)*3 - 1
            results.append((frame, occ, start_nt, end_nt))
    return results

    # Uso de map_all_proteins  "data/SARS-COV-2-MN908947.3.txt" y "data/seq-proteins.txt"
def map_all_proteins(genome_file, proteins_file):
    _, genome = read_fasta(genome_file)
    proteins = read_proteins(proteins_file)
    mapping = {}
    for pid, pseq in proteins.items():
        mapping[pid] = map_protein_to_genome(pseq, genome)
    return mapping

if __name__ == "__main__":
    import json, sys
    genome_file = sys.argv[1] if len(sys.argv)>1 else "data/SARS-COV-2-MN908947.3.txt"
    proteins_file = sys.argv[2] if len(sys.argv)>2 else "data/seq-proteins.txt"
    m = map_all_proteins(genome_file, proteins_file)
    print(json.dumps({k: v for k,v in m.items()}, indent=2))

# ------------------------------------------------------------------------

# FUNCIONES DEL PUNTO 4


# ------------------------------------------------------------------------


if __name__ == "__main__":
    print("===Analisis Genomico de SARS-COV-2===\n")

    print("Cargando archivos...\n")
    print("-"*50)
    files_data = load_all_files()

    print("🎯 PUNTO 1: Búsqueda de Genes")
    print("=" * 50)

    #gene_test()
    find_genes_index(files_data['genome']['MN908947.3'], files_data['genes'])

    print("\n🔄 PUNTO 2: Palíndromos | (Regiones Propensas a mutación)")
    print("=" * 50)

    #palindrome_test()
    analyze_palindromes_in_genes(files_data['genes'])

    print("\n🧬 PUNTO 3: MAPEO DE PROTEÍNAS")
    print("=" * 50)

    mapping = map_all_proteins("data/SARS-COV-2-MN908947.3.txt", "data/seq-proteins.txt")
    for prot, pos in mapping.items():
        print(f"{prot}: {pos if pos else 'NO encontrado'}")

    print("\n⚖️ PUNTO 4: COMPARACIÓN WUHAN vs TEXAS")
    print("=" * 50)