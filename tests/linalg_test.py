from src import linalg

def test_vectorize():
    mod5 = linalg.vectorize(lambda x: x % 5)
    assert mod5([4, 8, 2, 7, 11, 23]) == [4, 3, 2, 2, 1, 3]
    assert mod5([[11, 17], [-4, 10]]) == [[1, 2], [1, 0]]

def test_vector_mod():
    v = [6, 1, 4, 2]
    linalg.vector_mod(v, 3)
    assert linalg.vector_mod(v, 3) == [0, 1, 1, 2]

def test_add_vectors_mod():
    u = [10, 5, 7, 13, 19]
    v = [-3, -9, 4, 0, 12]
    assert linalg.sum_vectors(u, v) == [7, -4, 11, 13, 31]

def test_find_pivot():
    A = [[0, 3, 7],
         [0, 0, 1],
         [1, 0, 3],
         [0, 1, 7]]
    assert linalg.find_pivot(A, 0) == 2
    assert linalg.find_pivot(A, 1) == -1
    assert linalg.find_pivot(A, 2) == 2

def test_find_pivot2():
    A = [[0, 0, 3, 0],
         [0, 0, 0, 1],
         [0, 1, 0, 2]]
    assert linalg.find_pivot(A, 1) == 2
    assert linalg.find_pivot(A, 3) == -1

def test_matrix_mod():
    A = [[2, 10],
         [4, 13]]
    linalg.matrix_mod(A, 7)
    assert A == [[2, 3],
                 [4, 6]]

def test_swap():
    v = [1, 2, 3, 4, 5]
    linalg.swap(v, 0, 3)
    assert v == [4, 2, 3, 1, 5]

def test_transpose():
    A = [[0, 3, 7, 4],
         [0, 0, 1, 2],
         [1, 0, 3, 11]]
    assert linalg.transpose(A) == [[0, 0, 1],
                                   [3, 0, 0],
                                   [7, 1, 3],
                                   [4, 2, 11]]

def test_matrix_prod():
    A = [[1, 2, 3],
         [4, 5, 6]]
    B = [[1, 2],
         [3, 4],
         [5, 6]]
    assert linalg.matrix_prod(A, B) == [[22, 28],
                                        [49, 64]]

def test_naive_vector_prod():
    u = [1, 2, 3]
    v = [9, 8, 7]
    assert linalg.naive_vector_prod(u, v) == [9, 16, 21]

def test_rref():
    A = [[5, 2, 3],
         [2, 4, 1],
         [1, 0, 1]]
    assert linalg.rref(A) == [[5, 2, 3], [0.0, 3.2, -0.20000000000000018], [0.0, 0.0, 0.3749999999999999]]

def test_rref_assimetric1():
    A = [[5, 2],
         [2, 4],
         [1, 0]]
    assert linalg.rref(A) == [[5, 2], [0.0, 3.2], [0.0, 0.0]]

def test_rref_assimetric2():
    A = [[5, 2, 3],
         [2, 4, 1]]
    assert linalg.rref(A) == [[5, 2, 3], [0.0, 3.2, -0.20000000000000018]]

def test_echelon_mod_2():
    A = [[7, 3, 2],
         [3, 9, 1]]
    b = [3, 1]
    A, b = linalg.echelon_mod_2(A, b)
    assert A == [[1, 1, 0],
                 [0, 0, 1]]
    assert b == [1, 0]

def test_solve_mod_2():
    A = [[7, 3, 2],
         [3, 9, 1],
         [1, 6, 9]]
    b = [3, 1, 10]
    assert linalg.solve_mod_2(A, b) == [1, 0, 0]
