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
    palindrome_position = text.find(longest_palindrome)
    return {
        'sequence': longest_palindrome,
        'length': max_len,
        'position': palindrome_position,
    }

def palindrom_test():
    resultado = manacherAlgorithm("aaaa")
    print(f"Palíndromo más largo: {resultado}")
    
    # Probar con otros ejemplos
    print(f"Palíndromo más largo: {manacherAlgorithm("babad")}")  # Debería devolver "bab" o "aba"
    print(f"Palíndromo más largo: {manacherAlgorithm("cbbd")}")   # Debería devolver "bb"