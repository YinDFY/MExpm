import random
import gmpy2
from Crypto.Util.number import getPrime
gmpy2.get_context().precision = 2048




def calculate_average(numbers):
    if not numbers:
        return 0  # avoid divide zero

    total = sum(numbers)
    average = total / len(numbers)
    return average

def RandN(g1,p ,M ,flag = False,table_beta = [],table_alpha = [],n = 64,h = 64):
    """
    :param p: a large prime number
    :param M: a large number
    :param n: integer
    :param h: integer  h>=1
    :return: pair(x,X), where g^x = X mod p
    """
    times = gmpy2.mpz(gmpy2.div(p - 1, M))
    g = gmpy2.powmod(g1, times, p)
    #transform data into big number library form
    p = gmpy2.mpz(p)
    M = gmpy2.mpz(M)
    k = gmpy2.mpz(random.randint(1,n))
    h = gmpy2.mpz(h)
    n = gmpy2.mpz(n)
    #times = gmpy2.mpz(gmpy2.div(p-1,M))
    if flag == False:
        for i in range(n):
            alpha = gmpy2.mpz(random.randint(0,M-1))
            beta = gmpy2.powmod(g, alpha, p)
            table_alpha.append(alpha)
            table_beta.append(beta)
        return table_alpha,table_beta
    else:
        #Pair Generation
        x = gmpy2.mpz(0)
        X = gmpy2.mpz(1)
        i = 0
        while i < k:
           s = random.randint(0,n-1)
           xj = gmpy2.mpz(random.randint(1,h-1))
           x = gmpy2.add(x,gmpy2.mul(table_alpha[s],xj))
           x = gmpy2.f_mod(x,M)

           X = gmpy2.mul(X,gmpy2.powmod(table_beta[s],xj,p))
           X = gmpy2.f_mod(X, p)
           i = i + 1
           if gmpy2.f_mod(x,M) == 0:
                x = gmpy2.mpz(0)
                X = gmpy2.mpz(1)
                i = 0

        x = gmpy2.f_mod(x,M)
        X = gmpy2.f_mod(X,p)
        g = gmpy2.mpz(g)
        x = gmpy2.mpz(x)
        X = gmpy2.mpz(X)
        return g,x,X


# generate u
def generate_u(p,number):
    roots = set()  # Use a set to avoid duplicates
    while len(roots) < number:
        u = random.randint(2, p - 1)
        if gmpy2.gcd(u, p) == 1:
            roots.add(u)

    return list(roots)

# generate a
def generate_a(q, number):
    roots = set()  # Use a set to avoid duplicates
    while len(roots) < number:
        u = random.randint(2, q - 1)
        if gmpy2.gcd(u, q) == 1:
            roots.add(u)

    return list(roots)



def generate_large_prime(bit):
    prime_candidate = getPrime(bit)
    return gmpy2.mpz(prime_candidate)

def generate_p_q(bit):
    q = generate_large_prime(bit)  # produce bit-bit prime q
    i = gmpy2.mpz(2)
    while True:
        p_candidate = gmpy2.mul(q, i) + 1  # p = 2q + 1
        i = i + 1
        if gmpy2.is_prime(p_candidate):
            return p_candidate, q


if __name__ == "__main__":

    print(getPrime(32))
    print(getPrime(64))
    print(getPrime(128))
    print(getPrime(256))
    print(getPrime(256))
