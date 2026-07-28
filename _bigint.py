def inline_division(N: int, D: int) -> tuple[int, int]:
    # return integer result + remainder
    res = ""
    is_divisible = 0
    if D == 0:
        raise ValueError("Cannot divide by zero.")
    if N < D:
        return 0
    if N == D:
        return 1
    np = 0
    for dig in str(N):
        np = np * 10 + int(dig)
        qp = np // D
        res += str(qp)
        np = np % D
    return (int(res), np)

def big_exp_mod_N_aux(base: int, exp: int, mod: int) -> int:
    if exp == 1:
        return base % mod
    if exp == 0:
        return 1
    if exp % 2 == 0:
        half = big_exp_mod_N_aux(base, exp // 2, mod)
        return (half * half) % mod
    else:
        return (base * big_exp_mod_N_aux(base, exp - 1, mod)) % mod

def big_exp_mod_N(A: int, B: int, N: int) -> int:
    if N == 0:
        raise ValueError("Cannot perform mod operation with divisor 0.")
    if N == 1:
        return 0
    if A == 0 and B == 0:
        raise ValueError("0^0 is undefined.")
    if A == 0:
        return 0
    if B == 0:
        return 1
    return big_exp_mod_N_aux(A, B, N)

if __name__ == '__main__':
    while True:
        try:
            N = int(input("Insert the numerator: (0 to quit) "))
            if N == 0:
                break
            D = int(input("Insert the denominator: (0 to quit) "))
            if D == 0:
                break
            print(inline_division(N, D))
        except ValueError:
            print("Insert a valid number")