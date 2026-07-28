from sortedcontainers import SortedDict

from _bigint import inline_division

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

prime_map = sieve_of_eratosthenes(1000000)
known_primes = list[int]()
number_dict = SortedDict[int, int]()
node_arr = list[NumberGraphNode]()
n_nodes = 0
for i in range(1, len(prime_map)):
    if prime_map[i]:
        known_primes.append(i + 1)
        print(i + 1, end = " ")
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
        result = inline_division(num, prime)
        if result[1]:
            # if num is divisible by prime, add prime to list and pick the result
            res.append(prime)
            num = result[0]
        else:
            # num is not divisible by prime -> none of its divisors will be divisible by prime
            min_ind += 1
    if num in known_primes:
        res.append(num)
    return res

if __name__ == '__main__':
    num = -1
    try:
        num = input("Insert a number: ")
        num = int(num)
        for prime in find_known_factors(num):
            print(prime, end = " ")
        print()
    except ValueError as e:
        print(f"Error: {e.args[0]}")