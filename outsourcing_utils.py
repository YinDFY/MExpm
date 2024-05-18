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

def find_generator(L, phi_L):
    while True:
        g = gmpy2.mpz_random(gmpy2.random_state(), L - 1) + 1  # 1 <= g < L
        if gmpy2.gcd(g, L) != 1:
            continue  # Ensure g is coprime with L
        if gmpy2.powmod(g, phi_L, L) == 1:
            return g

def RandN(g1,p ,M ,flag = False,table_beta = [],table_alpha = [],n = 64,h = 64):
    """
    :param p: p = L*N 是两个大素数的乘法
    :param M: 是(L - 1)*(N - 1)是p的欧拉函数
    :param n: integer
    :param h: integer  h>=1
    :return: pair(x,X), where g^x = X mod p 要求X是模p阶为M的乘法群元素
    """
    #times = gmpy2.divexact(p - 1, M)
    g = g1

    #transform data into big number library form
    p = gmpy2.mpz(p)
    M = gmpy2.mpz(M)
    k = gmpy2.mpz(1)
    h = gmpy2.mpz(h)
    n = gmpy2.mpz(n)
    #times = gmpy2.mpz(gmpy2.div(p-1,M))
    if flag == False:
        for i in range(n):
            alpha = generate_u(M,1)[0]
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
           xj = generate_a(M,1)[0]
           xt = gmpy2.add(x, gmpy2.mul(table_alpha[s], xj))
           x = xt
           x = gmpy2.f_mod(x,M)
           X = gmpy2.mul(X,gmpy2.powmod(table_beta[s],xj,p))
           X = gmpy2.f_mod(X, p)
           i = i + 1
           if gmpy2.f_mod(x,M) == 0:
                x = gmpy2.mpz(0)
                X = gmpy2.mpz(1)
                i = 0

        x = generate_a(M,1)[0]
        X = gmpy2.powmod(g,x,p)

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
