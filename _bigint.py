def inline_division(N, D):
    res = ""
    last_index = len(D)
    if int(N[:last_index]) > int(D):
        last_index += 1
    np = int(N[:last_index])
    D = int(D)
    while True:
        res += str(np // D)
        np = np % D
        if last_index >= len(N):
            break
        np = np * 10 + int(N[last_index])
        last_index += 1
    return res

if __name__ == '__main__':
    N = -1
    D = -1
    while True:
        try:
            N = int(input("Insert the numerator: (0 to quit) "))
            if N == 0:
                break
            D = int(input("Insert the denominator: (0 to quit) "))
            if D == 0:
                break
            N = str(N)
            D = str(D)
            print(inline_division(N, D))
        except TypeError:
            print("Insert a valid number")