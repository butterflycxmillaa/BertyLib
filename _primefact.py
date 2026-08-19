import math
from sortedcontainers import SortedDict
from _bigint import big_exp_mod_N

class NumberGraphNode:
    def __init__(self, value: int):
        self.value = value

def sieve_of_eratosthenes(n: int) -> list[bool]:
    # Returns a list is_prime where is_prime[i] shows whether (i + 1) <is prime
    is_prime = [True] * n
    for num in range(2, int(n ** 0.5) + 1):
        if is_prime[num - 1]:
            # iterate over all the multiples of num
            for mul in range(num * 2, n + 1, num):
                is_prime[mul - 1] = False
    return is_prime

def segmented_soe(kp: list[int], L: int) -> list[int]:
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

prime_map = sieve_of_eratosthenes(1_000_000)
known_primes = list[int]()
number_dict = SortedDict[int, int]()
node_arr = list[NumberGraphNode]()
n_nodes = 0
for i in range(1, len(prime_map)):
    if prime_map[i]:
        known_primes.append(i + 1)
        print(i + 1, end = " ")
more_primes = segmented_soe(known_primes, 11_000_000)
for p in more_primes:
    print(p, end = " ")
print()

def is_num_present(num: int) -> bool:
    return num in number_dict

def return_num_index(num: int) -> int:
    return number_dict[num] if is_num_present(num) else -1

# returns the index of the node representing the number
# if it doesn't exist, it creates a new node and returns its index
def insert_num(num: int) -> int:
    if is_num_present(num):
        return number_dict[num]
    else:
        try:
            node_arr.append(NumberGraphNode(num))
            number_dict[num] = n_nodes
            n_nodes += 1
            return n_nodes - 1
        except ValueError as e:
            print(f"Invalid number: {num}. {e.args[0]}")
            return -1
  
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

if __name__ == '__main__':
    num = -1
    try:
        num = input("Insert a number: ")
        num = int(num)
        print(miller_rabin_primality(num))
    except ValueError as e:
        print(f"Error: {e.args[0]}")