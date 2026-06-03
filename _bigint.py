def inline_division(N: str, D: str) -> str:
    res = ""
    N = N.lstrip('0') or '0'
    D = D.lstrip('0') or '0'
    int_N = int(N)
    int_D = int(D)
    if D == '0':
        return "INF"
    if int_N < int_D:
        return "0"
    if int_N == int_D:
        return "1"
    np = 0
    for dig in N:
        np = np * 10 + int(dig)
        qp = np // int_D
        res += str(qp)
        np = np % int_D
    return res.lstrip('0') or '0'

if __name__ == '__main__':
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
        except ValueError:
            print("Insert a valid number")