from scipy.sparse import lil_matrix

def sparse_absolute(mat):
    """because numpy.absolute() forces it to be a dense matrix..."""
    mat = lil_matrix(mat)  # copy matrix
    # make each element absolute
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            mat[i, j] = abs(mat[i, j])
    return mat
