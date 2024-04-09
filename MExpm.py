from gmpy2 import *
from outsourcing_utils import *
import random
from Crypto.Util.number import getPrime

get_context().precision = 2048

def MExpm(u_list,a_list,p,N,table_alpha,table_beta,g1):
    blinding_pairs = []
    w = []
    w_ = []
    Aux = []
    L = mpz(mul(p, N))

    for i in range(4):
        g, x, X = RandN(g1, L, mul(p - 1, N - 1), True, table_beta, table_alpha)
        blinding_pairs.append([x, X])

    v1_inverse = invert(blinding_pairs[0][1], L)
    v2_inverse = invert(blinding_pairs[1][1], L)
    gtask = mpz(g)
    gveriy = mpz(g)
    a_sum = mpz(0)

    r = mpz(random.randint(2, N))
    Aux.append(f_mod(mul(r, blinding_pairs[2][1]), L))
    Aux.append(blinding_pairs[3][1])

    i = 0
    moudular_elur = mpz(mul(p - 1, N - 1))
    for u in u_list:
        wx = f_mod(mul(u, v1_inverse), L)
        w.append(wx)
        gtask = f_mod(mul(gtask, wx), L)
        a_sum = f_mod(add(a_sum, a_list[i]), moudular_elur)

        y = u
        wy = f_mod(mul(y, v2_inverse), L)
        w_.append(wy)
        gveriy = f_mod(mul(gveriy, wy), L)
        i = i + 1

    k1a_sum = f_mod(mul(blinding_pairs[0][0], a_sum), moudular_elur)
    t1z1 = f_mod(sub(k1a_sum, blinding_pairs[2][0]), moudular_elur)

    k2a_sum = f_mod(mul(blinding_pairs[1][0], a_sum), moudular_elur)
    t2z2 = f_mod(sub(k2a_sum, blinding_pairs[3][0]), moudular_elur)

    m = []
    m_ = []
    qu = generate_a(moudular_elur, 2)
    for a in a_list:
        mx = f_mod(sub(a, t1z1), moudular_elur)
        m.append(f_mod(mul(mx, invert(qu[0], moudular_elur)), moudular_elur))

        my = f_mod(sub(a, t2z2), moudular_elur)
        m_.append(f_mod(mul(my, invert(qu[1], moudular_elur)), moudular_elur))

    ans_x = []
    ans_y = []

    TK = blinding_pairs[2][1]
    VK = blinding_pairs[3][1]

    mid_TK = powmod(gtask, t1z1, L)

    mid_VK = powmod(gveriy, t2z2, L)

    for j in range(len(w)):
        ans_x.append(powmod(w[j], m[j], L))
        ans_y.append(powmod(w_[j], m_[j], L))

    prod_x = mpz(1)
    prod_y = mpz(1)
    for j in range(len(w)):
        prod_x = f_mod(mul(prod_x, ans_x[j]), L)
        prod_y = f_mod(mul(prod_y, ans_y[j]), L)
    prod_x = powmod(prod_x, qu[0], L)
    prod_y = powmod(prod_y, qu[1], L)

    TK = f_mod(mul(TK, prod_x), L)
    VK = f_mod(mul(VK, prod_y), L)

    TK = f_mod(mul(TK, mid_TK), L)
    VK = f_mod(mul(VK, mid_VK), L)

    return TK, VK, r

if __name__ == "__main__":
    bitlong1 = 512
    bitlong2 = 1024
    number = 1001 #batch size

    #MExpm
    print("|-----------------------------------MExpm-----------------------------------|")
    p = getPrime(bitlong1)
    N = getPrime(bitlong1)
    L = gmpy2.mpz(gmpy2.mul(p, N))  #L is  bitlong2
    u1 = generate_u(N, number)
    g1 = mpz(random.randint(2, L))
    table_alpha, table_beta = RandN(g1, L, mul(p - 1, N - 1), False)

    k = mpz(random.randint(2, N))
    kN = mpz(mul(k, N))

    y = []

    for u in u1:
        y1 = f_mod(add(u, kN), L)
        y.append(y1)
    a = generate_a(gmpy2.mul(p - 1, N - 1), number)
    TK, VK, r = MExpm(y, a, p, N, table_alpha, table_beta, g1)

    #compute real answer
    real_answer = mpz(1)
    for i in range(number):
        temp = powmod(u1[i], a[i], N)
        real_answer = f_mod(mul(real_answer, temp), N)

    recovery_answer = f_mod(TK, N)
    if TK == VK:
        print("Verification is passed")
    if recovery_answer == real_answer:
        print("Real answer is true")



    #Exp  Exp operations in numbers
    print("|-----------------------------------Exp-----------------------------------|")
    p, q = generate_p_q(bitlong2)
    u1 = generate_u(p, number)
    a1 = generate_a(q, number)
    u = []
    t = mpz(1)
    for tmp, tee in zip(u1, a1):
        qe = powmod(tmp, tee, p)
        t = f_mod(mul(t, qe), p)
        u.append(powmod(tmp, tee, p))
