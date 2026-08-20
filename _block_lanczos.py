import random
 
def popcount(x):
    return bin(x).count('1')
 
def mat_vec_gf2(A_rows, v, n):
    r = 0
    for i in range(n):
        if popcount(A_rows[i] & v) & 1:
            r |= (1 << i)
    return r
 
def block_mat_vec(A_rows, V, n):
    return [mat_vec_gf2(A_rows, v, n) for v in V]
 
def block_inner(V, W):
    b = len(V)
    res = [0] * b
    for i in range(b):
        row = 0
        for j in range(b):
            if popcount(V[i] & W[j]) & 1:
                row |= (1 << j)
        res[i] = row
    return res
 
def gf2_matmul_square(X, Y, b):
    res = [0] * b
    for i in range(b):
        r = 0
        for k in range(b):
            if (X[i] >> k) & 1:
                r ^= Y[k]
        res[i] = r
    return res
 
def block_right_mult(V, C, b):
    result = [0] * b
    for j in range(b):
        acc = 0
        for k in range(b):
            if (C[k] >> j) & 1:
                acc ^= V[k]
        result[j] = acc
    return result
 
def invert_gf2_square(matrix, b):
    aug = [matrix[i] | (1 << (b + i)) for i in range(b)]
    for c in range(b):
        piv = next((i for i in range(c, b) if (aug[i] >> c) & 1), None)
        if piv is None:
            return None
        aug[c], aug[piv] = aug[piv], aug[c]
        for i in range(b):
            if i != c and (aug[i] >> c) & 1:
                aug[i] ^= aug[c]
    mask = (1 << b) - 1
    return [(aug[i] >> b) & mask for i in range(b)]
 
def gf2_rref_kernel(rows, ncols):
    nrows = len(rows)
    pivots = []
    r = 0
    for c in range(ncols):
        piv = next((i for i in range(r, nrows) if (rows[i] >> c) & 1), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(nrows):
            if i != r and (rows[i] >> c) & 1:
                rows[i] ^= rows[r]
        pivots.append(c)
        r += 1
        if r == nrows:
            break
    pivot_set = set(pivots)
    basis = []
    for free in (c for c in range(ncols) if c not in pivot_set):
        vec = 1 << free
        for i, pc in enumerate(pivots):
            if (rows[i] >> free) & 1:
                vec |= (1 << pc)
        basis.append(vec)
    return basis
 
def reduce_to_basis(vectors):
    rows = [v for v in vectors if v != 0]
    r = 0
    n_used = 0
    while r < len(rows):
        # trova colonna pivot per la riga r
        col = rows[r].bit_length() - 1
        if col < 0:
            rows.pop(r)
            continue
        for i in range(len(rows)):
            if i != r and (rows[i] >> col) & 1:
                rows[i] ^= rows[r]
        r += 1
    return [v for v in rows if v != 0]
 
def block_lanczos_kernel(A_rows, n, block_size=8, seed=None):
    b = min(block_size, n)
    rnd = random.Random(seed)
    max_iters = 10 * (n // b + 5)
    V_prev = [0] * b
    V_curr = [rnd.getrandbits(n) for _ in range(b)]
    Tinv_prev = None
    krylov_vectors = []
    deflations = []
    expansions = 0
    for it in range(1, max_iters + 1):
        if all(v == 0 for v in V_curr):
            break
        U = block_mat_vec(A_rows, V_curr, n)
        T = block_inner(V_curr, U)
        Tinv = invert_gf2_square(T, b)
        if Tinv is None:
            for z in gf2_rref_kernel(T, b):
                cand = 0
                for k in range(b):
                    if (z >> k) & 1:
                        cand ^= V_curr[k]
                if cand != 0 and mat_vec_gf2(A_rows, cand, n) == 0:
                    deflations.append(cand)
            V_curr = [rnd.getrandbits(n) for _ in range(b)]
            V_prev = [0] * b
            Tinv_prev = None
            continue
        krylov_vectors.extend(V_curr)
        expansions += 1
        W = block_mat_vec(A_rows, U, n)
        X = block_inner(V_curr, W)
        C = gf2_matmul_square(Tinv, X, b)
        if Tinv_prev is not None:
            Y = block_inner(V_prev, W)
            D = gf2_matmul_square(Tinv_prev, Y, b)
            VprevD = block_right_mult(V_prev, D, b)
        else:
            VprevD = [0] * b
        VcurrC = block_right_mult(V_curr, C, b)
        V_next = [U[k] ^ VcurrC[k] ^ VprevD[k] for k in range(b)]
        V_prev, V_curr = V_curr, V_next
        Tinv_prev = Tinv
    Av_rows = [mat_vec_gf2(A_rows, v, n) for v in krylov_vectors]
    m = len(Av_rows)
    track = [1 << i for i in range(m)]
    M = Av_rows[:]
    pr = 0
    for col in range(n):
        piv = next((r for r in range(pr, m) if (M[r] >> col) & 1), None)
        if piv is None:
            continue
        M[piv], M[pr] = M[pr], M[piv]
        track[piv], track[pr] = track[pr], track[piv]
        for r in range(m):
            if r != pr and (M[r] >> col) & 1:
                M[r] ^= M[pr]
                track[r] ^= track[pr]
        pr += 1
        if pr == m:
            break
    kernel_vectors = list(deflations)
    for r in range(m):
        if M[r] == 0:
            v = 0
            for k in range(m):
                if (track[r] >> k) & 1:
                    v ^= krylov_vectors[k]
            if v != 0 and mat_vec_gf2(A_rows, v, n) == 0:
                kernel_vectors.append(v)
    return reduce_to_basis(kernel_vectors)
 
def find_row_dependencies(M_rows, block_size=8, seed=None):
    m = len(M_rows)
    A_rows = [0] * m
    for i in range(m):
        row = 0
        for j in range(m):
            if popcount(M_rows[i] & M_rows[j]) & 1:
                row |= (1 << j)
        A_rows[i] = row
    candidates = block_lanczos_kernel(A_rows, m, block_size=block_size, seed=seed)
    dependencies = []
    for v in candidates:
        combo = [i for i in range(m) if (v >> i) & 1]
        xor_check = 0
        for i in combo:
            xor_check ^= M_rows[i]
        if xor_check == 0:
            dependencies.append(combo)
    return dependencies