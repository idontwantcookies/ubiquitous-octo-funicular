from collections.abc import Iterable
from numbers import Number
from typing import Callable

import sympy


Vector = list[Number]
Matrix = list[Vector]

def vectorize(f: Callable[[Number], Number]) -> Callable[[Iterable | Number], Iterable | Number]:
    def g(u: Iterable | Number) -> Iterable | Number:
        if isinstance(u, Iterable):
            acc = []
            for x in u:
                acc.append(g(x))
            return acc
        return f(u)
    return g

def vector_mod(v: Vector, p: int) -> Vector:
    u = []
    for x in v:
        u.append(x % p)
    return u

def matrix_mod(A: Matrix, p: int):
    A = A.copy()
    for i, row in enumerate(A):
        A[i] = vector_mod(row, p)
    return A

def swap(A: list, i: int, j: int):
    A[i], A[j] = A[j], A[i]

def find_pivot(A: Matrix, j: int) -> int:
    i, n, m = j, len(A), len(A[0])
    while i < n and i < m and A[i][j] == 0:
        i += 1
    if i in (m, n): return -1
    return i

def sum_vectors(*vectors: list[Vector]) -> Vector:
    acc = []
    for row in zip(*vectors):
        acc.append(sum(row))
    return acc

def scale_vector(u: Vector, alpha: Number) -> Vector:
    w = []
    for x in u:
        w.append(x * alpha)
    return w

def naive_vector_prod(u: Vector, v: Vector) -> Vector:
    w = []
    for x, y in zip(u, v):
        w.append(x * y)
    return w

def transpose(A: Matrix) -> Matrix:
    N = len(A)
    M = len(A[0])
    T = [[0 for _ in range(N)] for _ in range(M)]
    for i in range(N):
        for j in range(M):
            T[j][i] = A[i][j]
    return T

def matrix_prod(A: Matrix, B: Matrix):
    prod = []
    B = transpose(B)
    N, M = len(A), len(B)
    for i in range(N):
        prod.append([])
        for j in range(M):
            prod[i].append(sum(naive_vector_prod(A[i], B[j])))
    return prod

def gauss_reduce_row(row: Vector, pivot: Vector, i: int):
    m = - 1 / pivot[i] * row[i]
    row = sum_vectors(row, scale_vector(pivot, m))
    return row

def rref(A: Matrix) -> Matrix:
    '''Reduz uma matriz N x M à sua forma escalonada (Reduced Row Echelon Form) usando o método
    de Gauss. Complexidade: O(N²)'''
    A = A.copy()
    N = len(A)
    for i in range(N):
        p = find_pivot(A, i)
        if p == -1 or A[p][p] == 0: continue
        swap(A, i, p)
        pivot_row = A[i]
        for j in range(i + 1, N):
            A[j] = gauss_reduce_row(A[j], pivot_row, i)
    return A

def kernel(A: sympy.Matrix | Matrix) -> Matrix:
    A = sympy.Matrix(A)
    ker: list[sympy.Matrix] = A.nullspace()
    return [list(u.transpose()) for u in ker]
