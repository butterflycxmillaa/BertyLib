from sortedcontainers import SortedDict

from _bigint import inline_division

class NumberGraphNode:
    def __init__(self, value: str):
        if (len(value) > 1 and not value[1:].isdigit()) or not value[0].isdigit():
            raise ValueError("Value must be a string representation of a number.")
        if int(value) <= 1:
            raise ValueError("Value must be greater than 1.")
        self.value = value

class PrimeFactorization:
    n_nodes = 0
    number_dict = SortedDict[str, int]()
    node_arr = list[NumberGraphNode]()
    known_primes = list[str]()

    def __init__(self):
        known_primes = self.sieve_of_eratosthenes(1000000)
        for i in range(len(known_primes)):
            if known_primes[i]:
                self.known_primes.append(str(i))

    @staticmethod
    def sieve_of_eratosthenes(n: int) -> list[bool]:
        is_prime = [True] * (n + 1)
        is_prime[0] = is_prime[1] = False
        for i in range(2, int(n**0.5) + 1):
            if is_prime[i]:
                for j in range(i * 2, n + 1, i):
                    is_prime[j] = False
        return is_prime

    def is_num_present(cls, num: str) -> bool:
        return num in cls.number_dict

    # returns the index of the node representing the number
    # if it doesn't exist, it creates a new node and returns its index
    def insert_num(cls, num: str) -> int:
        if cls.is_num_present(num):
            return cls.number_dict[num]
        else:
            try:
                cls.node_arr.append(NumberGraphNode(num))
                cls.number_dict[num] = cls.n_nodes
                cls.n_nodes += 1
                return cls.n_nodes - 1
            except ValueError as e:
                print(f"Invalid number: {num}. {e.args[0]}")
                return -1

if __name__ == '__main__':
    pf = PrimeFactorization()
    count = 0
    for prime in pf.known_primes:
        count += 1
    print(count)