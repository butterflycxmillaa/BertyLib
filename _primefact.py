import math
from sortedcontainers import SortedDict
from _bigint import big_exp_mod_N

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

kp = set_upper_lim_kp(10_000_000)
for p in kp:
    print(p, end = " ")
print()

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
            exp /= 2
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

def pollard_rho(num: int) -> int:
    # won't work if num is prime
    # instead it will enter an infinite loop (fix needed)
    if num % 2 == 0:
        return 2
    while True:
        c = 1
        x0 = 2
        def apply_polinomial(x: int):
            return (x ** 2 + c) % num
        # sets tortoise and hare pointers
        T = H = x0
        while True:
            # execute the iteration
            T = apply_polinomial(T)
            H = apply_polinomial(apply_polinomial(H))
            gcd = euclidean_gcd_algorithm(abs(T - H), num)
            if 1 < gcd < num:
                return gcd
            if gcd == num:
                break
        c += 1

def compute_qs_params(num: int) -> tuple[int, int]:
    # computes the length of the optimal factor base
    ln_num = math.log(num)
    ln_ln_num = math.log(ln_num)
    c = 1 / math.isqrt(2)
    ln_B = c * math.sqrt(ln_num * ln_ln_num)
    B = math.exp(ln_B)
    k = int((B / math.log(B)) / 2)
    digits = len(str(num))
    if digits <= 30:
        optimal_k = int(300 + (digits / 30) ** 2 * 200)
    elif digits <= 60:
        optimal_k = int(500 + ((digits - 30) / 30) ** 2.2 * 8500)
    elif digits <= 90:
        optimal_k = int(9000 + ((digits - 60) / 30) ** 2.4 * 250000)
    else:
        optimal_k = int(260000 + ((digits - 90) / 10) ** 2.5 * 740000)
    if digits >= 40:
        optimal_k = int(k * 0.65)
    optimal_B = int(2 * optimal_k * math.log(max(optimal_k, 2)))
    # the first k is the ideal one
    # the second k is the one used in practice
    # same goes for B
    return (k, optimal_k, B, optimal_B)

def legendre_symbol(num: int, p: int) -> bool:
    return pow(num, (p - 1) // 2, p) == 1

if __name__ == '__main__':
    num = -1
    try:
        num = input("Insert a number: ")
        num = int(num)
        print(miller_rabin_primality(num))
    except ValueError as e:
        print(f"Error: {e.args[0]}")