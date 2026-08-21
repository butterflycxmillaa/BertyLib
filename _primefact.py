import math
from _bigint import big_exp_mod_N
from _block_lanczos import find_row_dependencies
from sortedcontainers import SortedDict

from _TimeLimitExceeded_exception import TimeLimitExceededError as tlee

known_primes = []
upper_lim = 0

def sieve_of_eratosthenes(n: int) -> list[int]:
    # returns the list of primes p such that p <= n
    is_prime = [True] * n
    for num in range(2, int(n ** 0.5) + 1):
        if is_prime[num - 1]:
            # iterate over all the multiples of num
            for mul in range(num * 2, n + 1, num):
                is_prime[mul - 1] = False
    return [i + 1 for i in range(len(is_prime)) if is_prime[i] and i != 0]

def segmented_soe(kp: list[int], L: int) -> list[int]:
    if not kp:
        return sieve_of_eratosthenes(L)
    M = kp[-1] + 1
    if L < M: return []
    # candidates[0] points to M
    # candidates[-1] points to L
    dim = L - M + 1
    candidates = [True] * dim
    max_p = math.isqrt(L)
    kp_segm = [p for p in kp if p <= max_p]
    for p in kp_segm:
        start = max((M + p - 1) // p * p, p ** 2)
        for i in range(start - M, dim, p):
            candidates[i] = False
    return [M + i for i in range(dim) if candidates[i]]

def set_upper_lim_kp(lim: int) -> list[int]:
    global known_primes, upper_lim
    if upper_lim == 0:
        sqrt_lim = int(math.isqrt(lim) + 1)
        known_primes = sieve_of_eratosthenes(sqrt_lim)
    expand_primes = segmented_soe(known_primes, lim)
    known_primes += expand_primes
    upper_lim = lim
    return known_primes

kp = set_upper_lim_kp(1_000_000)

def find_known_factors(num: int) -> list[int]:
    # factorizes using the known primes and returns a list of the factors found
    res = list[int]()
    # defines the already explored known primes
    min_ind = 0
    while num not in known_primes and num != 1 and min_ind < len(known_primes):
        # check if it's divisible by the first available prime
        prime = known_primes[min_ind]
        result = (num // prime, num % prime)
        if result[1] == 0:
            # if num is divisible by prime, add prime to list and pick the result
            res.append(prime)
            num = result[0]
        else:
            # num is not divisible by prime -> none of its divisors will be divisible by prime
            min_ind += 1
    if num in known_primes:
        res.append(num)
    return res

def miller_rabin_primality(num: int) -> bool:
    exp = num - 1
    A = 100
    final = big_exp_mod_N(A, exp, num)
    if final == 1:
        while exp % 2 == 0:
            exp //= 2
            # perform the big exp calculation once again
            new = big_exp_mod_N(A, exp, num)
            if final == 1:
                if new == 1 or new == num - 1:
                    final = new
                    continue
                return False
        return True
    return False

def euclidean_gcd_algorithm(A: int, B: int) -> int:
    A = abs(A)
    B = abs(B)
    # swap A and B so that A >= B
    if A == 0:
        return B
    if B == 0:
        return A
    if B > A:
        A, B = B, A
    result = (A // B, A % B)
    while result[1] != 0:
        A = B
        B = result[1]
        result = (A // B, A % B)
    return B

def pollard_rho(num: int, max_c: int = 10, max_iters: int = 10 ** 6) -> int:
    # won't work if num is prime
    # instead it will enter an infinite loop (fix needed)
    if num % 2 == 0:
        return 2
    c = 1
    while c <= max_c:
        x0 = 2
        def apply_polinomial(x: int):
            return (x ** 2 + c) % num
        # sets tortoise and hare pointers
        iters = 0
        T = H = x0
        while iters <= max_iters:
            # execute the iteration
            T = apply_polinomial(T)
            H = apply_polinomial(apply_polinomial(H))
            gcd = euclidean_gcd_algorithm(abs(T - H), num)
            if 1 < gcd < num:
                return gcd
            if gcd == num:
                break
            iters += 1
        c += 1
    raise tlee("No factor was found within set iteration limit")

def generate_N_primes(N: int) -> list[int]:
    primes = [2]
    if N <= 0: return []
    if N == 1: return primes
    idx = 0
    num = 3
    while idx < N - 1:
        if miller_rabin_primality(num):
            primes.append(num)
            idx += 1
        num += 2
    return primes

def expand_N_primes(kp: list[int], N: int) -> list[int]:
    len_kp = len(kp)
    if N <= len_kp:
        return kp[:N]
    last_prime = kp[-1] + 2
    n_primes = N - len_kp
    idx = 0
    while idx < n_primes:
        if miller_rabin_primality(last_prime):
            kp.append(last_prime)
            idx += 1
        last_prime += 2
    return kp

def compute_qs_params(num: int) -> tuple[int, int]:
    # computes the length of the optimal factor base
    digits = len(str(num))
    if digits < 10:
        optimal_k = 2 + 1.2 * (digits ** 1.4)
    elif digits <= 30:
        optimal_k = 35 + ((digits - 10) / 20) ** 2 * 265
    elif digits <= 60:
        optimal_k = 300 + ((digits - 30) / 30) ** 2.2 * 8700
    elif digits <= 90:
        optimal_k = 9000 + ((digits - 60) / 30) ** 2.4 * 251000
    else:
        optimal_k = 260000 + ((digits - 90) / 10) ** 2.5 * 740000
    optimal_k = int(optimal_k)
    if digits >= 40:
        optimal_k = int(optimal_k * 0.65)
    optimal_k = max(optimal_k, 5)
    optimal_B = int(2 * optimal_k * math.log(max(2 * optimal_k, 2)))
    return (optimal_k, optimal_B)

def legendre_symbol(num: int, p: int) -> bool:
    return p == 2 or pow(num, (p - 1) // 2, p) == 1

def generate_factor_base(num: int, k: int, B: int):
    global known_primes, upper_lim

    if upper_lim < B:
        set_upper_lim_kp(B)
    factor_base = []
    for p in known_primes:
        if len(factor_base) - 1 >= k:
            break
        if legendre_symbol(num, p):
            factor_base.append(p)
    current_lim = upper_lim
    while len(factor_base) - 1 < k:
        current_lim *= 2
        set_upper_lim_kp(current_lim)
        for p in known_primes:
            if p <= known_primes[len(factor_base) - 1]:
                continue
            if len(factor_base) - 1 >= k:
                break
            if legendre_symbol(num, p):
                factor_base.append(p)
    return factor_base

def perform_sieving(qx: int, fb: int) -> list[int]:
    # tries to divide qx by all the factors inside of the factor base
    # if qx is not fb-smooth, then returns list [-1] * k
    k = len(fb)
    res = [0] * k
    for factor_n, factor in enumerate(fb):
        while qx % factor == 0:
            res[factor_n] += 1
            qx //= factor
    if qx == 1: return res
    else: return [-1] * k

def quadratic_sieve(num: int) -> int:
    k, B = compute_qs_params(num)
    fb = generate_factor_base(num, k, B)
    len_fb = len(fb)
    x = math.isqrt(num) + 1
    x_arr = [-1] * (len_fb + 1)
    qx_arr = [-1] * (len_fb + 1)
    fact = [0] * (len_fb + 1)
    idx = 0
    while idx < len_fb + 1:
        qx = (x ** 2 - num)
        sieve = perform_sieving(qx, fb)
        if sieve[0] == -1:
            # the Q(x) is not a B-smooth value, therefore x should be discarded
            x += 1
            continue
        # Q(x) here is guaranteed to be B-smooth, so it should be added to the matrix
        x_arr[idx] = x
        qx_arr[idx] = qx
        exp_arr_mod_2 = [exp % 2 for exp in sieve]
        for pos, exp in enumerate(exp_arr_mod_2[::-1]):
            fact[idx] += (1 << pos) * exp
        x += 1
        idx += 1
    deps = find_row_dependencies(fact, 64)
    for dep in deps:
        # compute values of X and Y
        X = 1
        Y = 1
        for elem in dep:
            X *= x_arr[elem]
            Y *= qx_arr[elem]
        X %= num
        Y = math.sqrt(Y)
        g = euclidean_gcd_algorithm(X - Y, num)
        if g == 1 or abs(g) == num: continue
        return int(g)

def factorize_num_aux(num: int, factors: list[int]) -> None:
    global upper_lim

    if miller_rabin_primality(num):
        factors.append(num)
        return
    if num < upper_lim:
        kfs = find_known_factors(num)
        factors += kfs
        return
    g = 1
    try:
        g = pollard_rho(num)
    except tlee:
        g = quadratic_sieve(num) 
    finally:
        # g is a factor of num, so factorize both g and num // g
        factorize_num_aux(g, factors)
        factorize_num_aux(num // g, factors)
    return

def factorize_num(num: int) -> SortedDict[int, int]:
    res = SortedDict[int, int]()
    if num == 0:
        res[0] = 1
        return res
    else:
        res[0] = 0
    if num < 0:
        res[-1] = 1
        num *= -1
    elif num > 0:
        res[-1] = 0
    if num == 1:
        return res
    factors = []
    factorize_num_aux(num, factors)
    for f in factors:
        if f not in res.keys():
            res[f] = 0
        res[f] += 1
    return res

if __name__ == '__main__':
    num = -1
    try:
        num = 465782945632
        print(f"{num} =", end = " ")
        res = factorize_num(num)
        for idx, fact in enumerate(res.keys()):
            if fact > 0:
                print(f"({fact} ^ {res[fact]}){' *' if fact != res.keys()[-1] else ''}", end = " ")
        print()
    except ValueError as e:
        print(f"Error: {e.args[0]}")