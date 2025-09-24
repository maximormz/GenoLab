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