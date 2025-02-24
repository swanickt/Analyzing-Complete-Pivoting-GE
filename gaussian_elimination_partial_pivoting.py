import copy

def partial_piv_gaussian_elimination(A: list[list[float]], b: list[float]) -> list[float]:
    """
    This function implements Gaussian Elimination with partial pivoting (i.e. only
    row pivoting). In particular, we use the largest entry on or below the diagonal in
    each column as the pivot. A pivot vector rpiv is used to track the row ordering,
    without having to actually swap them in the matrix.

    Preconditions:
        - A is a square matrix.
        - If A is nxn, then b is an nx1 vector.
    """
    A = copy.deepcopy(A)
    b = copy.deepcopy(b)
    n = len(A)
    x = [0.0] * n # the vector we will return, the solution to the system Ax = b.
    rpiv = list(range(n)) # A row pivot vector so we can use indirect addressing

    ### Gaussian Elimination ###
    for i in range(n-1): # loop over the diagonal of the matrix
        # find the row index of the element with the largest absolute value in the (i)th column of A.
        max_row = max(range(i, n), key=lambda l: abs(A[rpiv[l]][i]))
        if A[rpiv[max_row]][i] == 0:
            print("A is a singular matrix.")
            return None

        # update the pivot vector (simulates swapping rows i and max_row)
        rpiv[i], rpiv[max_row] = rpiv[max_row], rpiv[i]

        for j in range(i+1, n): # compute multipliers below the current pivot
            multiplier = A[rpiv[j]][i] / A[rpiv[i]][i]
            for k in range(i, n): # R_j <- (R_j - multiplier * R_i)
                A[rpiv[j]][k] -= (multiplier * A[rpiv[i]][k])
            b[rpiv[j]] -= (multiplier * b[rpiv[i]]) # update the RHS of Ax = b

    # We now have an upper triangular system

    ### Backwards Substitution to solve for x ###
    for i in range(n-1, -1, -1): # loop backwards over columns, from n-1 to 0 (inclusive).
        if A[rpiv[i]][i] == 0:
            print("Singular matrix encountered during back substitution! :(")
            return None
        summation = 0
        for j in range(i+1, n):
            summation += x[j] * A[rpiv[i]][j]
        x[i] = (b[rpiv[i]] - summation) / A[rpiv[i]][i]

    return x

if __name__ == '__main__':
    # 3x3 toy example with a non-singular matrix A:
    A_test = [[ 2, 4, -2],
              [ 4, 14, 0],
              [-1, 10, 7]]

    b_test = [2, -1, 1]

    x = partial_piv_gaussian_elimination(A_test, b_test)
    print("Using the partial pivoting GE function:")
    print(f" {A_test[0]}[x_1]   [ {b_test[0]}]          "
          f"[x_1]   [{x[0]:.3f}]\n {A_test[1]}[x_2] = [{b_test[1]}]  yields  [x_2] = [  {x[1]:.3f}]"
          f"\n{A_test[2]}[x_3]   [ {b_test[2]}]          "
          f"[x_3]   [ {x[2]:.3f}]\n")