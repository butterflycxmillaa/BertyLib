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
        for i in range(2, int(n ** 0.5) + 1):
            if is_prime[i]:
                for j in range(i * 2, n + 1, i):
                    is_prime[j] = False
        return is_prime

    def is_num_present(cls, num: str) -> bool:
        return num in cls.number_dict
    
    def return_num_index(cls, num: str) -> int:
        return cls.number_dict[num] if cls.is_num_present(num) else -1

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

    @classmethod       
    def find_known_factors(cls, num: str) -> list[str]:
        res = list[str]()
        return cls.find_known_factors_aux(num, res, 0)
    
    @classmethod
    def find_known_factors_aux(cls, num: str, res: list[str], threshold_min_ind: int) -> list[str]:
        if num == "1": return list[str]
        if num in cls.known_primes:
            res.append(num)
            return res
        threshold_max = int(int(num) ** 0.5) + 1
        for known_prime in cls.known_primes[threshold_min_ind::]:
            if int(known_prime) <= threshold_max:
                result = inline_division(num, known_prime)
                print(f"{num} / {known_prime} = {result[0]}, full division: {result[1]}")
                if result[1]:
                    res.append(known_prime)
                    return cls.find_known_factors_aux(result[0], res, threshold_min_ind)
                else:
                    threshold_min_ind += 1
            else: break
        if num != "1":
            res.append(num)
        return res

if __name__ == '__main__':
    pf = PrimeFactorization()
    num = -1
    try:
        num = input("Insert a number: ")
        if len(num) > 7:
            print("Number is too long")
        else:
            num = int(num)
            num = str(num)
    except ValueError:
        print("Number is not valid")
    for prime in pf.find_known_factors(num):
        print(prime, end = " ")
    print()