Vector = list[int]
Matrix = list[Vector]

def vector_mod(v: Vector, p: int):
    u = []
    for x in v:
        u.append(x % p)
    return u

def swap(A: list, i: int, j: int):
    A[i], A[j] = A[j], A[i]

def find_pivot(A: Matrix, j: int):
    i, n = j, len(A)
    while i < n and A[i][j] == 0:
        i += 1
    return i

def sum_vectors(u: Vector, v: Vector):
    if len(u) != len(v): raise ValueError("Vectors' length must be the same.")
    w = []
    for x, y in zip(u, v):
        w.append(x + y)
    return w

def naive_vector_prod(u: Vector, v: Vector):
    w = []
    for x, y in zip(u, v):
        w.append(x * y)
    return w

def matrix_mod(A: Matrix, p: int):
    for i, row in enumerate(A):
        A[i] = vector_mod(row, p)

def transpose(A: Matrix):
    N = len(A)
    M = len(A[0])
    T = [[0 for _ in range(N)] for _ in range(M)]
    for i in range(N):
        for j in range(M):
            T[j][i] = A[i][j]
    return T

def matrix_prod(A, B):
    prod = []
    B = transpose(B)
    N, M = len(A), len(B)
    for i in range(N):
        prod.append([])
        for j in range(M):
            prod[i].append(sum(naive_vector_prod(A[i], B[j])))
    return prod

def echelon_mod_2(A: Matrix, b: Vector):
    A = A.copy()
    matrix_mod(A, 2)
    b = vector_mod(b, 2)
    n = len(A)
    for i in range(n):
        p = find_pivot(A, i)
        if p >= n: continue
        if A[p][p] == 0: continue
        swap(A, i, p)
        swap(b, i, p)
        for j in range(i + 1, n):
            if j > len(A): break
            if A[j][i] == 0: continue
            A[j][i:] = sum_vectors(A[i][i:], A[j][i:])
            A[j][i:] = vector_mod(A[j][i:], 2)
            b[j] = (b[i] + b[j]) % 2
    return A, b

def solve_mod_2(A: Matrix, b: Vector):
    A, b = echelon_mod_2(A, b)
    A = transpose(A)
    A, b = echelon_mod_2(A, b)
    return b
