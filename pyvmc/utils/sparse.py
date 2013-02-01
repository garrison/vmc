from scipy.sparse import lil_matrix

def sparse_absolute(mat):
    """because numpy.absolute() forces it to be a dense matrix..."""
    mat = lil_matrix(mat)  # copy matrix
    # make each element absolute
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            mat[i, j] = abs(mat[i, j])
    return mat

def is_hermitian(mat):
    # For dense matrices you can simply test (mat == mat.H), but for some
    # reason this always returns False for sparse matrices.
    if mat.shape[0] != mat.shape[1]:
        return False
    for i in range(mat.shape[0]):
        for j in range(i, mat.shape[1]):
            if mat[(i, j)] != mat[(j, i)].conjugate():
                return False
    return True
