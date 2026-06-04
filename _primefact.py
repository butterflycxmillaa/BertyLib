from sortedcontainers import SortedDict

class NumberGraphNode:
    def __init__(self, value: str):
        self.value = value

class PrimeFactorization:
    n_nodes = 0
    number_dict = SortedDict[str, int]()
    node_arr = list[NumberGraphNode]()

    def is_num_present(cls, num: str) -> bool:
        return num in cls.number_dict
    
    def insert_num(cls, num: str) -> int:
        if cls.is_num_present(num):
            return cls.number_dict[num]
        else:
            cls.node_arr.append(NumberGraphNode(num))
            cls.number_dict[num] = cls.n_nodes
            cls.n_nodes += 1
            return cls.n_nodes - 1

if __name__ == '__main__':
    pf = PrimeFactorization()
    print(pf.insert_num('2'))
    print(pf.insert_num('3'))
    print(pf.insert_num('2'))
    print(pf.number_dict)