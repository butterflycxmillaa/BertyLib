def float_to_fraction(x: str, sep: str = '.') -> str:
    try:
        float(x)
    except ValueError:
        return "Invalid input: not a valid floating-point number"
    if sep not in x:
        return x + '/1'

    integer_part, fractional_part = x.split(sep)
    denominator = 10 ** len(fractional_part)
    try:
        integer_part = int(integer_part)
        fractional_part = int(fractional_part)
        numerator = int(integer_part) * denominator + int(fractional_part)
    except ValueError:
        return "Invalid input: not a valid floating-point number"
    
    return f'{numerator}/{denominator}'

if __name__ == '__main__':
    while True:
        try:
            x = input("Insert a floating-point number (or '0' to quit): ")
            if x == '0':
                break
            print(float_to_fraction(x))
        except ValueError:
            print("Insert a valid floating-point number")