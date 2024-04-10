from gmpy2 import *
from outsourcing_utils import *
import random
from Crypto.Util.number import getPrime
import time
get_context().precision = 2048


def MExpm(u_list,a_list,L,q,N,t,table_alpha,table_beta,g1):
    blinding_pairs = []
    w = []
    w_ = []
    # keyGen
    start_time = time.time()
    for i in range(4):
        g, x, X = RandN(g1, L, mul(p - 1, N - 1), True, table_beta, table_alpha)
        blinding_pairs.append([x, X])

    v1_inverse = invert(blinding_pairs[0][1], L)  # 1mi
    v2_inverse = invert(blinding_pairs[1][1], L)  # 1mi

    gtask = mpz(g)
    gveriy = mpz(g)
    a_sum = mpz(0)
    elur = mpz(mul(p - 1, N - 1))  # 1mm
    r = mpz(random.randint(2, N))
    qu = generate_a(elur, 2)
    i = 0

    m_index = random.randint(1, len(u_list) - 1)
    d = invert(a_list[m_index], elur)  # 1mi
    rd = powmod(r, d, L)  # 1exp

    for u in u_list:
        wx = f_mod(mul(u, v1_inverse), L)  # nmm
        w.append(wx)
        gtask = f_mod(mul(gtask, wx), L)  # nmm
        a_sum = f_mod(add(a_sum, a_list[i]), elur)
        if i == m_index:
            u = f_mod(mul(u, rd), L)  # 1mm
        wy = f_mod(mul(u, v2_inverse), L)  # nmm
        w_.append(wy)
        gveriy = f_mod(mul(gveriy, wy), L)  # nmm
        i = i + 1

    k1a_sum = f_mod(mul(blinding_pairs[0][0], a_sum), elur)  # 1mm
    t1z1 = f_mod(sub(k1a_sum, blinding_pairs[2][0]), elur)
    t1 = qu[0]
    z1 = f_mod(mul(t1z1, invert(t1, elur)), elur)  # 1mm+1minv

    k2a_sum = f_mod(mul(blinding_pairs[1][0], a_sum), elur)  # 1mm
    t2z2 = f_mod(sub(k2a_sum, blinding_pairs[3][0]), elur)
    t2 = qu[0]
    z2 = f_mod(mul(t2z2, invert(t2, elur)), elur)  # 1mm+1minv

    m = []
    m_ = []
    for a in a_list:
        mx = f_mod(sub(a, t1z1), elur)
        m.append(mx)
        # m.append(mx)

        my = f_mod(sub(a, t2z2), elur)
        m_.append(my)
        # m_.append(my)
    TK = f_mod(mul(blinding_pairs[2][1], r), L)  # 1mm
    end_time = time.time()
    keygen_time = (end_time - start_time) * 1000
    """
    Compute
    """
    start_time = time.time()

    mid_TK = powmod(gtask, z1, L)
    mid_VK = powmod(gveriy, z2, L)

    ans_x = []
    ans_y = []
    for j in range(len(w)):
        ans_x.append(powmod(w[j], m[j], L))
        ans_y.append(powmod(w_[j], m_[j], L))
    end_time = time.time()
    compute_time = (end_time - start_time) * 1000
    """
    Verify
    """
    start_time = time.time()
    VK = blinding_pairs[3][1]
    prod_x = mpz(1)
    prod_y = mpz(1)
    for j in range(len(w)):
        prod_x = f_mod(mul(prod_x, ans_x[j]), L)
        prod_y = f_mod(mul(prod_y, ans_y[j]), L)
    mid_TK = powmod(mid_TK, t1, L)
    mid_VK = powmod(mid_VK, t2, L)

    prod_x = f_mod(mul(mid_TK, prod_x), L)
    prod_x = f_mod(mul(TK, prod_x), L)

    prod_y = f_mod(mul(mid_VK, prod_y), L)
    prod_y = f_mod(mul(VK, prod_y), L)
    end_time = time.time()
    verify_time = (end_time - start_time) * 1000

    return prod_x, prod_y, r#, keygen_time, compute_time, verify_time
if __name__ == "__main__":

    bitlong1 = 512
    bitlong2 = 1024
    number = 1000
    start_time = time.time()
    p = getPrime(bitlong1)
    N = getPrime(bitlong1)
    L = mul(p, N)
    u1 = generate_u(L, number)
    a1 = generate_a(gmpy2.mul(p - 1, N - 1), number)
    g1 = mpz(random.randint(2, L - 1))
    table_alpha, table_beta = RandN(g1, L, mul(p - 1, N - 1), False)

    end_time2 = time.time()
    k = mpz(random.randint(2, p))
    kN = mpz(mul(k, N))

    y = []
    a = []
    for u in u1:
        y1 = f_mod(add(u, kN), L)
        y.append(y1)

    for tmp in a1:
        a.append(tmp)
    end_time = time.time()
    TK, VK, r = MExpm(y, a, L, p, N, N, table_alpha, table_beta, g1)

    if TK == VK:
        print("Verification is passed")

    start_time = time.perf_counter()
    r_inverse = invert(r, L)  # 1minv
    TK = f_mod(mul(TK, r_inverse), L)  # 1mm

    #compute real answer using direct exp operations
    real_answer = mpz(1)
    for base,exponent in zip(u1,a):
        tmp = powmod(base,exponent,N)
        real_answer = f_mod(mul(real_answer,tmp),N)

    if real_answer == f_mod(TK,N):
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
