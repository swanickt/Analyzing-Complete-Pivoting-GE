import copy

def no_piv_gaussian_elimination(A: list[list[float]], b: list[float]) -> list[float]:
    """
    This function implements Gaussian Elimination with no pivoting. Thus, it terminates
    when a zero pivot is encountered on the diagonal. Since we are not performing
    any row or column interchanges in this version, we do not need any pivot vectors.

    Preconditions:
        - A is a square matrix.
        - If A is nxn, then b is an nx1 vector.
    """
    A = copy.deepcopy(A) # Avoid mutating the inputs
    b = copy.deepcopy(b)
    n = len(A)
    x = [0.0] * n # the vector we will return, the solution to the system Ax = b.

    ### Gaussian Elimination ###
    for i in range(n-1): # loop over the diagonal of the matrix
        pivot = A[i][i]
        if pivot == 0:
            print("Zero pivot encountered in no pivot function! :(")
            return None
        for j in range(i+1, n): # compute multipliers below the current pivot
            multiplier = A[j][i] / pivot
            for k in range(i, n): # R_j <- (R_j - multiplier * R_i)
                A[j][k] -= (multiplier * A[i][k])
            b[j] -= (multiplier * b[i]) # update the RHS of Ax = b

    # We now have an upper triangular system

    ### Backwards Substitution to solve for x ###
    for i in range(n-1, -1, -1): # loop backwards over columns, from n-1 to 0 (inclusive).
        if A[i][i] == 0: # stop if matrix is singular
            print("Singular matrix encountered during back substitution! :(")
            return None
        summation = 0
        for j in range(i+1, n):
            summation += x[j] * A[i][j]
        x[i] = (b[i] - summation) / A[i][i]

    return x

if __name__ == '__main__':
    # 3x3 toy example with a non-singular matrix A:
    A_test = [[ 2, 4, -2],
              [ 4, 14, 0],
              [-1, 10, 7]]

    b_test = [2, -1, 1]


    x = no_piv_gaussian_elimination(A_test, b_test)
    print("Using the no pivoting GE function:")
    print(f" {A_test[0]}[x_1]   [ {b_test[0]}]          "
          f"[x_1]   [{x[0]:.3f}]\n {A_test[1]}[x_2] = [{b_test[1]}]  yields  [x_2] = [  {x[1]:.3f}]"
          f"\n{A_test[2]}[x_3]   [ {b_test[2]}]          "
          f"[x_3]   [ {x[2]:.3f}]\n")