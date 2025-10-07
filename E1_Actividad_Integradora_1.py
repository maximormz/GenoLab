"""
COMENTARIOS DEL CÓDIGO SARS-COV-2

Este código realiza análisis genómico del virus SARS-COV-2, incluyendo:
1. Búsqueda de genes en genomas
2. Detección de palíndromos en secuencias genéticas
3. Mapeo de proteínas en el genoma
4. Comparación entre diferentes cepas del virus

INSTRUCCIÓN: No modificar la lógica, solo agregar comentarios
"""

import os
import re

# CONSTANTES
# Obtener el directorio actual del script
current_dir = os.path.dirname(os.path.abspath(__file__))

# Tabla de codones para traducción de ADN a aminoácidos
# Cada codon de 3 nucleótidos se mapea a un aminoácido específico
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
    # Construir ruta completa al archivo
    file_path = os.path.join(current_dir,filename)
    normalized_path = os.path.normpath(file_path)

    genome = ''

    try:
        # Abrir y leer el archivo
        with open(normalized_path, 'r') as file:
            for line in file:
                # Saltar líneas de cabecera (empiezan con '>')
                if not line.startswith('>'):
                    # Concatenar secuencia eliminando saltos de línea
                    genome += line.strip("\n")

        return genome
    except FileNotFoundError:
        print(f"Error: El archivo '{filename}' no se encontró en {normalized_path}.")
        return ""

# Cargar archivos de genes
def load_gene_file(filename):
    # Similar a load_genome_file pero para archivos de genes
    file_path = os.path.join(current_dir,filename)
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
    file_path = os.path.join(current_dir,filename)
    normalized_path = os.path.normpath(file_path)

    # Diccionario para almacenar proteínas: {nombre: secuencia}
    proteins = {}
    protein_name = None

    try:
        with open(normalized_path, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    # Nueva proteína - obtener nombre (eliminar '>' y espacios)
                    protein_name = line[1:].strip()
                    proteins[protein_name] = ''
                elif protein_name is not None and not line.startswith('>'):
                    # Concatenar secuencia de aminoácidos
                    proteins[protein_name] += line.strip()

        return proteins
    except FileNotFoundError:
        print(f"Error: El archivo '{filename}' no se encontró en {normalized_path}.")
        return ""

# Retornar todos los archivos cargados
def load_all_files():
    # Cargar dos genomas de SARS-COV-2
    genome_data_1 = load_genome_file('SARS-COV-2-MN908947.3.txt')
    genome_data_2 = load_genome_file('SARS-COV-2-MT106054.1.txt')
    
    # Cargar secuencias de genes específicos
    gene_M_data = load_gene_file('gen-M.txt')
    gene_S_data = load_gene_file('gen-S.txt')
    gene_ORF1AB_data = load_gene_file('gen-ORF1AB.txt')
    
    # Cargar secuencias de proteínas
    protein_data = load_protein_file('seq-proteins.txt')

    # Retornar estructura organizada de datos
    return {
        'genome': {
            'Wuhan_2019' : genome_data_1,
            'Texas_2020' : genome_data_2},
        'genes': {
            'M': gene_M_data,
            'S': gene_S_data,
            'ORF1AB': gene_ORF1AB_data
        },
        'proteins': protein_data
    }

# ------------------------------------------------------------------------

# FUNCIONES DE APOYO
def kmp_search_all_occurrences(text, pattern): # Retorna un generador (Convertir en lista al recibirlo)
    # Algoritmo Knuth-Morris-Pratt para búsqueda eficiente de patrones
    if not text or not pattern:
        return
    
    n = len(text)
    m = len(pattern)

    i,j = 0,0

    # Precomputar tabla de fallos
    failure_table = build_failure_table(pattern)

    # Búsqueda principal
    while i < n:
        if text[i] == pattern[j]:
            i += 1
            j += 1

            if j == m:
                # Patrón encontrado - yield posición de inicio
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
    # LPS = Longest Prefix Suffix
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
    # Transformar texto para algoritmo Manacher
    # Inserta caracteres especiales entre cada carácter
    n = len(text)
    e = 2 * n + 3  # Longitud de la nueva cadena

    new_string = ["$"] * e
    new_string[0] = "@"   # Carácter inicial
    new_string[-1] = "#"  # Carácter final

    for i in range(n):
        new_string[2 * i + 2] = text[i]

    return e, new_string

def manacherAlgorithm(text):
    # Algoritmo para encontrar el palíndromo más largo
    centro = limite = 0
    e,string = new_phrase(text)
    # Array para almacenar longitudes de palíndromos
    P = [0] * e

    for i in range(1,e - 1):
        if i < limite:
            simetrica = 2 * centro - i
            P[i] = min(limite - i,P[simetrica])

        # Expandir alrededor del centro i
        gap = P[i] + 1
        while string[i-gap] == string[i+gap]:
            P[i] += 1
            gap += 1
        
        # Actualizar centro y límite si es necesario
        if i + P[i] > limite:
            limite = i + P[i]
            centro = i

    # Encontrar el palíndromo más largo
    max_len = max(P)
    center_index = P.index(max_len)

    # Calcular posición en el texto original
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
    # Buscar genes dentro del genoma usando KMP
    if not genome_sequence:
        print("⚠️  Advertencia: La secuencia del genoma está vacía")
        return {}
    if not gene_sequence:
        print("⚠️  Advertencia: No hay secuencias de genes para buscar")
        return {}
    
    results = {}

    # Buscar cada gen en el genoma
    for gene_name,sequence in gene_sequence.items():
        if not sequence:
            print(f"⚠️  Advertencia: Secuencia vacía para el gen {gene_name}")
            continue
        
        # Buscar todas las ocurrencias del gen
        direct_indexes = list(kmp_search_all_occurrences(genome_sequence, sequence))

        # Almacenar resultados
        results[gene_name] = {
            'indexes': direct_indexes,
            'sequence': sequence,
            'first_12': sequence[:12],  # Primeros 12 nucleótidos
            'last_12': sequence[-12:],  # Últimos 12 nucleótidos
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

# ------------------------------------------------------------------------

# FUNCIONES DEL PUNTO 2
def save_palindrome_to_file(gene_name, palindrome_sequence): 
    # Guardar resultados de palíndromos en archivos
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
    # Mostrar resultados de análisis de palíndromos
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
    # Analizar palíndromos en todas las secuencias de genes
    results = {}
    for gene_name, gene_sequence in gene_sequence.items():
        print(f"\n🔬 Analizando palíndromos en Gen {gene_name}...")

        # Encontrar palíndromo más largo
        longest_palindrome = manacherAlgorithm(gene_sequence)
        file_name = save_palindrome_to_file(gene_name,longest_palindrome)

        # Almacenar resultados
        results[gene_name] = {
            'gene_sequence': gene_sequence,
            'longest_palindrome': longest_palindrome['sequence'],
            'palindrome_length': longest_palindrome['length'],
            'palindrome_position': longest_palindrome['position'],
            'saved_to_file' : file_name
        }

        display_palindrome_results(gene_name,results[gene_name])

    return results

# ------------------------------------------------------------------------

# FUNCIONES DEL PUNTO 3
# Tabla de codones para traducción (similar a la anterior pero con formato diferente)
CODON_TABLE_mapper = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '', 'TAG': '',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def read_fasta(path):
    # Leer archivo FASTA - versión simplificada
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
    # Limpiar secuencia - solo mantener ATCGN
    seq = re.sub('[^ATCGN]', '', seq)
    return header, seq

def read_proteins(path):
    # Leer archivo de proteínas en formato FASTA
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
                    # Guardar proteína anterior
                    seq = ''.join(cur_lines).upper()
                    seq = re.sub('[^A-Z]', '', seq)
                    proteins[cur_id] = seq
                # Nueva proteína
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
    # Traducir ADN a proteína dado un marco de lectura
    aa = []
    for i in range(frame, len(dna) - 2, 3):
        codon = dna[i:i+3]
        # Obtener aminoácido o 'X' si no se encuentra
        aa.append(CODON_TABLE_mapper.get(codon, 'X'))
    return ''.join(aa)

def kmp_search_all(text, pattern):
    # Implementación KMP para búsqueda en texto
    if not pattern:
        return []
    m = len(pattern)
    lps = [0] * m
    length = 0
    i = 1
    # Construir tabla LPS
    while i < m:
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

    res = []
    n = len(text)
    i = j = 0
    # Búsqueda principal
    while i < n:
        if text[i] == pattern[j]:
            i += 1
            j += 1
            if j == m:
                res.append(i - j)
                j = lps[j - 1]
        else:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1
    return res

def verify_slippery(genome, nt_start, aa_target, max_aa=120):
    # Verificar sitios "slippery" - posibles frameshifts
    i_aa = 0
    shift = False
    while i_aa < min(max_aa, len(aa_target)):
        codon = genome[nt_start:nt_start + 3]
        if len(codon) < 3:
            break
        aa = CODON_TABLE_mapper.get(codon, "X")
        if aa == aa_target[i_aa]:
            # Coincidencia - avanzar normal
            nt_start += 3
            i_aa += 1
            continue
        if not shift:
            # Primer desfase - posible frameshift
            nt_start += 1
            shift = True
            continue
        break
    return i_aa

def map_protein_to_genome(protein_seq, genome_seq):
    # Mapear secuencia proteica al genoma
    protein_seq = protein_seq.upper()
    results = []
    
    # Buscar en los 3 marcos de lectura
    for frame in range(3):
        translated = translate_dna(genome_seq, frame)
        occs = kmp_search_all(translated, protein_seq)
        for occ in occs:
            # Calcular posiciones en nucleótidos
            start_nt = frame + occ * 3 + 1
            end_nt = start_nt + len(protein_seq) * 3 - 1
            results.append((frame, occ, start_nt, end_nt))

    # Si no se encuentra, buscar con mecanismo slippery
    if not results:
        pref = protein_seq[:8]
        for frame in range(3):
            i = frame
            while i <= len(genome_seq) - 3:
                match_len = verify_slippery(genome_seq, i, pref, max_aa=120)
                if match_len >= len(pref) - 1:
                    start_nt = i + 1
                    end_nt = start_nt + match_len * 3 - 1
                    occ = i // 3
                    results.append((frame, occ, start_nt, end_nt))
                    return results  # Retornar el primer match con slippery
                else:
                    i += 3
    return results

def print_information(mapping):
    # Imprimir resultados del mapeo de proteínas
    total_found = sum(1 for v in mapping.values() if v)
    total_notfound = sum(1 for v in mapping.values() if not v)
    print(f"Total proteinas: {len(mapping)} | Encontradas: {total_found} | No encontradas: {total_notfound}")
    print("-------------------------------------------------------")
    print(f"{'🧬 ID':<8} {' 🧪 Frame':<8} {' 🏁 Inicio':<9} {' 🏁 Fin':<10} {'📏 Longitud':<12} {'✅ Estado'}")
    print("-" * 65)

    for pid, results in mapping.items():
        if results:
            for frame, occ, start_nt, end_nt in results:
                length = end_nt - start_nt + 1
                print(f"{pid:<15} {frame:<7} {start_nt:<10} {end_nt:<10} {length:<10} {'OK'}")
        else:
            print(f"{pid:<15} {'-':<7} {'-':<10} {'-':<10} {'-':<10} {'NO'}")

    print("-" * 65)
    print("=======================================================")

def map_all_proteins(WuhanGenome, proteins_file):
    # Mapear todas las proteínas al genoma de Wuhan
    proteins = read_proteins(proteins_file)
    mapping = {}
    for pid, pseq in proteins.items():
        mapping[pid] = map_protein_to_genome(pseq, WuhanGenome)
    return print_information(mapping)

# ------------------------------------------------------------------------

# FUNCIONES DEL PUNTO 4
def compare_genomes_table(Genome_Name,WuhanGenome, TexasGenome, max_difs=50):
    """Compara dos genomas y muestra diferencias en tabla organizada"""
    seq1 = WuhanGenome
    seq2 = TexasGenome

    difs = []
    n = min(len(seq1), len(seq2))

    # Encontrar todas las diferencias
    for i in range(n):
        if seq1[i] != seq2[i]:
            # Encontrar el codón completo que contiene la diferencia
            codon_i = i - (i % 3)
            codon1 = seq1[codon_i:codon_i + 3]
            codon2 = seq2[codon_i:codon_i + 3]
            # Traducir a aminoácidos
            aa1 = CODON_TABLE.get(codon1, "X") if len(codon1) == 3 else "X"
            aa2 = CODON_TABLE.get(codon2, "X") if len(codon2) == 3 else "X"
            difs.append((i, codon1, aa1, codon2, aa2))

    # Mostrar resultados en tabla
    print(f"\nDiferencias entre {Genome_Name[0]} y {Genome_Name[1]} (max {max_difs}): {len(difs)}\n")
    print(f"{'Pos':>5}  {Genome_Name[0]:^10}  AA  {Genome_Name[1]:^10}  AA")
    print("-" * 40)

    for d in difs[:max_difs]:
        print(f"{d[0]:5d}  {d[1]:^10}  {d[2]}  {d[3]:^10}  {d[4]}")

    if len(difs) > max_difs:
        print(f"... y {len(difs) - max_difs} más diferencias")

    return difs

# ------------------------------------------------------------------------

if __name__ == "__main__":
    print("\n===Analisis Genomico de SARS-COV-2===\n")

    print("Cargando archivos...\n")
    print("-"*50)
    files_data = load_all_files()

    # PUNTO 1: Búsqueda de genes en el genoma
    print("🎯 PUNTO 1: Búsqueda de Genes")
    print("=" * 50)
    find_genes_index(files_data['genome']['Wuhan_2019'], files_data['genes'])
    
    # PUNTO 2: Análisis de palíndromos en genes
    print("\n🔄 PUNTO 2: Palíndromos | (Regiones Propensas a mutación)")
    print("=" * 50)
    analyze_palindromes_in_genes(files_data['genes'])

    # PUNTO 3: Mapeo de proteínas al genoma
    print("\n🧬 PUNTO 3: MAPEO DE PROTEÍNAS")
    print("=" * 50)
    map_all_proteins(files_data["genome"]['Wuhan_2019'], "seq-proteins.txt")

    # PUNTO 4: Comparación entre cepas
    print("\n⚖️  PUNTO 4: COMPARACIÓN WUHAN vs TEXAS")
    print("=" * 50)
    compare_genomes_table(list(files_data['genome'].keys()),files_data['genome']['Wuhan_2019'],files_data['genome']['Texas_2020'])
    print()