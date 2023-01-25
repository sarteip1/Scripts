gmpy2 = None
import math
import random

import gmpy2

# Seed the random number generator
random.seed(0)


# Define an exception class for when the modular inverse cannot be calculated
class InvError(Exception):
    def __init__(self, v):
        self.value = v


# Define a function to convert inputs to integers
def Int(x):
    # If gmpy2 is not available, return the standard integer representation of x
    if gmpy2 is None:
        return int(x)
    # Otherwise, return the mpz representation of x using gmpy2
    else:
        return gmpy2.mpz(x)


# Define a function to calculate the modular inverse of a mod n
def inv(a, n):
    # Reduce a modulo n
    a %= n
    # If gmpy2 is not available, use the built-in pow function to calculate the inverse
    if gmpy2 is None:
        try:
            return pow(a, -1, n)
        # If the inverse cannot be calculated, raise an exception
        except ValueError:
            import math
            raise InvError(math.gcd(a, n))
    # If gmpy2 is available, use the gcdext function to calculate the inverse
    else:
        g, s, t = gmpy2.gcdext(a, n)
        # If the gcd is not 1, raise an exception
        if g != 1:
            raise InvError(g)
        # Otherwise, return the modular inverse
        return s % n


# Define a class to represent points on an elliptic curve
class ECpoint(object):
    def __init__(self, A, B, N, x, y, *, prepare=True):
        # If prepare is True, perform some basic checks and reductions
        if prepare:
            N = Int(N)
            A, B, x, y = [Int(e) % N for e in [A, B, x, y]]
            if (y ** 2 - x ** 3 - A * x - B) % N != 0:
                raise ValueError
        # Assign the input values to the object
        self.A, self.B, self.N, self.x, self.y = A, B, N, x, y

    # Define the addition operation for the class
    def __add__(self, other):
        A, B, N = self.A, self.B, self.N
        Px, Py, Qx, Qy = self.x, self.y, other.x, other.y
        # If the two points being added are the same, use a different formula for the slope
        if Px == Qx and Py == Qy:
            s = ((Px * Px * 3 + A) * inv(Py * 2, N)) % N
        # Otherwise, use the standard formula for the slope
        else:
            s = ((Py - Qy) * inv(Px - Qx, N))
        x = (s * s - Px - Qx) % N
        y = (s * (Px - x) - Py) % N
        # Return a new ECpoint object
        return ECpoint(A, B, N, x, y, prepare=False)

    # Define the multiplication operation for the class
    def __rmul__(self, other):
        other = Int(other - 1)
        r = self
        while True:
            if other & 1:
                r = r + self
                if other == 1:
                    return r
            other >>= 1
            self = self + self


# Define a function to perform a binary search
def BinarySearch(f, a, b):
    while a < b:
        m = (a + b) // 2
        if f(m):
            b = m
        else:
            a = m + 1
    # Assert that the function is true at the final value
    assert a == b and f(a), (a, b, f(a))
    return a


# Define a function to calculate the minimum bit size of a factor that can be found with a given number of curves
def FactorBitSize(bound, ncurves, miss_prob=0.1):
    import math
    bound_log2 = math.log2(bound)

    # Define a function to calculate the probability of finding a factor of a given bit size
    def Prob(factor_log2):
        x = factor_log2 / bound_log2
        return x ** -x

    def F(factor_log2):
        return (1. - Prob(factor_log2)) ** ncurves >= miss_prob

    # Perform a binary search to find the minimum bit size of a factor that can be found with the given number of curves
    return BinarySearch(lambda x: F(x / 1000.), math.log2(bound) * 1000., 512 * 1000.) / 1000.


# Define a function to calculate the minimum number of curves needed to find a factor of a given bit size
def NeededCurves(bound, target_fac_log2, miss_prob=0.1):
    def F(ncurves):
        return FactorBitSize(bound, ncurves, miss_prob) >= target_fac_log2

    return round(BinarySearch(lambda x: F(x), 1, 10 ** 15))


def Work(bound, bound_pow, *, logs=[0.], cache={}):
    if bound not in cache:
        import math
        for ilog in range(100):
            if get_prime(1 << ilog) >= bound:
                break
        cnt_primes = BinarySearch(lambda i: get_prime(i) >= bound, 0, 1 << ilog)
        bound_pow = max(bound, bound_pow)
        bound_log = math.log2(bound)
        bound_pow_log = math.log2(bound_pow)
        while cnt_primes >= len(logs):
            plog = math.log2(get_prime(len(logs)))
            sum_plog = plog
            while True:
                sum_plog += plog
                if sum_plog >= max(bound_log, bound_pow_log):
                    sum_plog -= plog
                    break
            logs.append(logs[-1] + sum_plog)
        cache[bound] = logs[cnt_primes]
    return cache[bound]


def Work2(bound, bound_pow, factor_log):
    return Work(bound, bound_pow) * NeededCurves(bound, factor_log)


def OptimalBoundLog(factor_log, bound_pow):
    import math
    # Initialize the minimum work and bound values
    mwork = None
    bound_log_start = 10
    # Iterate through possible bound log values
    for bound_log in range(bound_log_start, math.floor(factor_log) + 1):
        bound = 2 ** bound_log
        # Calculate the work needed for this bound value
        work = Work(bound, bound_pow) * NeededCurves(bound, factor_log)
        # Check if the work is greater than the current minimum work, if so, break the loop
        if mwork is not None and work > mwork[0]:
            break
        # Update the minimum work and bound values if necessary
        if mwork is None or work < mwork[0]:
            mwork = (work, bound_log)
    else:
        bound_log = bound_log_start + 1
    # Fine-tune the bound log value
    mult = 200
    mwork = None
    for bound_log2 in range((bound_log - 3) * mult, bound_log * mult + 1):
        bound = round(2 ** (bound_log2 / mult))
        work = Work(bound, bound_pow) * NeededCurves(bound, factor_log)
        if mwork is None or work < mwork[0]:
            mwork = (work, bound_log2 / mult)
    return mwork[1]



def get_prime(i, *, primes=[2, 3]):
    while i >= len(primes):
        for n in range(primes[-1] + 2, 1 << 62, 2):
            isp = True
            for p in primes:
                if p * p > n:
                    break
                if n % p == 0:
                    isp = False
                    break
            if isp:
                primes.append(n)
                break
    return primes[i]


def prime_power(idx, bound, bound_pow, *, cache={}):
    key = (bound, bound_pow)
    if key not in cache:
        bound_pow = max(bound, bound_pow)
        r = []
        for i in range(1 << 62):
            p = get_prime(i)
            if p >= bound:
                break
            m = p
            while True:
                m2 = m * p
                if m2 >= bound_pow:
                    break
                m = m2
            r.append(m)
        cache[key] = r
    return cache[key][idx] if idx < len(cache[key]) else None


def ExcInfo(ex):
    return f'{type(ex).__name__}: {ex}'


def ProcessCurve(*, n, bound, bound_pow, shared, rand_seed, curve_idx, num_curves):
    try:
        random.seed(rand_seed)
        x0, y0, a = [random.randrange(1, n) for i in range(3)]
        b = (y0 ** 2 - x0 ** 3 - a * x0) % n

        P = ECpoint(a, b, n, x0, y0)
        for i in range(1 << 62):
            if shared.get('finish', False):
                return {'ex': 'Interrupted: Finishing...'}
            k = prime_power(i, bound, bound_pow)
            if i > 0 and i % 2000 == 0 or k is None:
                Print(('.', 'o')[k is None], end='', flush=True)
            if k is None:
                break
            P = k * P
    except InvError as e:
        d = e.value
        if d != n:
            return {'factors': sorted([d, n // d])}
    except BaseException as ex:
        shared['finish'] = True
        return {'ex': ExcInfo(ex)}
    return {}


def Print(*pargs, **nargs):
    print(*pargs, **{'flush': True, **nargs})


def ECM(n):
    print(f'Factoring {n}')

    if fermat_prime(n):
        return [n]

    bound = max(int(math.e ** (1 / 2 * math.sqrt(math.log(n) * math.log(math.log(n))))), 100)
    bound_pow = max(bound, 1 << 18)
    factor_log = math.log2(n) / 2
    max_curves = NeededCurves(bound, factor_log)

    try:
        ncurves, report_time = 0, 0
        shared = {}
        for icurve in range(1 << 62):
            result = ProcessCurve(n=n, bound=bound, bound_pow=bound_pow, shared=shared,
                                  rand_seed=random.randrange(1 << 48), curve_idx=icurve, num_curves=max_curves)
            if 'factors' in result:
                return result['factors']
            ncurves += 1
    except BaseException as ex:
        print(f'\nException: {ExcInfo(ex)}. Finishing, wait!')
    finally:
        shared['finish'] = True
        print()

    return [n]


def fermat_prime(n, trials=32):
    if n <= 16:
        return n in (2, 3, 5, 7, 11, 13)
    for i in range(trials):
        if pow(random.randint(2, n - 2), n - 1, n) != 1:
            return False
    return True


def gen_random_prime(bits):
    while True:
        n = random.randrange(1 << (bits - 1), 1 << bits)
        if fermat_prime(n):
            return n


def Prod(it):
    import functools
    return functools.reduce(lambda x, y: x * y, it, 1)


def run():
    pi = 1569275433846670190958947524742355825262149596324592390121
    fs = ECM(pi)
    Print('Factors:', ', '.join([('C', 'P')[fermat_prime(e)] + f' {e}' for e in fs]))


if __name__ == '__main__':
    run()
