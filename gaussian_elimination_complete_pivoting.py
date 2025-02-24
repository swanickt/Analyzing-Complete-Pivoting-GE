import copy

def complete_piv_gaussian_elimination(A: list[list[float]], b: list[float]) -> list[float]:
    """
    This function implements Gaussian Elimination with complete pivoting (i.e. both
    row and column pivoting to make each pivot on the diagonal as large as possible).
    We search the entire remaining unreduced sub-matrix for the largest entry and permute
    it into the diagonal pivot position. Pivot vectors rpiv and cpiv are used to track the row
    and column orderings, respectively. They facilitate indirect addressing of the modified
    matrix at each step of the GE process.

    Preconditions:
        - A is a square matrix.
        - If A is nxn, then b is an nx1 vector.
    """
    A = copy.deepcopy(A)
    b = copy.deepcopy(b)
    n = len(A)
    x = [0.0] * n # the vector we will return, the solution to the system Ax = b.
    rpiv = list(range(n)) # The row pivot vector to use indirect addressing
    cpiv = list(range(n)) # The column pivot vector to use indirect addressing

    ### Gaussian Elimination ###
    for i in range(n-1): # loop over the diagonal of the matrix

        # find the row and column index of the element with the largest absolute value
        # in the sub matrix with top left corner at (i, i).
        k, j = max(((k, j) for k in range(i, n) for j in range(i, n)), key=lambda kj: abs(A[rpiv[kj[0]]][cpiv[kj[1]]]))

        if A[rpiv[k]][cpiv[j]] == 0:
            print("A is a singular matrix.")
            return None
        rpiv[i], rpiv[k] = rpiv[k], rpiv[i] # update the row pivot vector
        cpiv[i], cpiv[j] = cpiv[j], cpiv[i]  # update the column pivot vector

        for m in range(i+1, n): # compute multipliers below the current pivot
            multiplier = A[rpiv[m]][cpiv[i]] / A[rpiv[i]][cpiv[i]]
            for l in range(i, n): # R_m <- (R_m - multiplier * R_i)
                A[rpiv[m]][cpiv[l]] -= (multiplier * A[rpiv[i]][cpiv[l]])
            b[rpiv[m]] -= (multiplier * b[rpiv[i]]) # update the RHS of Ax = b

    # We now have an upper triangular system

    ### Backwards Substitution to solve for x ###
    for p in range(n - 1, -1, -1):  # loop backwards over columns, from n-1 to 0 (inclusive).
        if A[rpiv[p]][cpiv[p]] == 0:
            print("Singular matrix encountered during back substitution! :(")
            return None
        summation = 0
        for j in range(p + 1, n):
            summation += x[cpiv[j]] * A[rpiv[p]][cpiv[j]]
        x[cpiv[p]] = (b[rpiv[p]] - summation) / A[rpiv[p]][cpiv[p]]

    return x

if __name__ == '__main__':
    # 3x3 toy example with a non-singular matrix A:
    A_test = [[ 2, 4, -2],
              [ 4, 14, 0],
              [-1, 10, 7]]

    b_test = [2, -1, 1]

    x = complete_piv_gaussian_elimination(A_test, b_test)
    print("Using the complete pivoting GE function:")
    print(f" {A_test[0]}[x_1]   [ {b_test[0]}]          "
          f"[x_1]   [{x[0]:.3f}]\n {A_test[1]}[x_2] = [{b_test[1]}]  yields  [x_2] = [  {x[1]:.3f}]"
          f"\n{A_test[2]}[x_3]   [ {b_test[2]}]          "
          f"[x_3]   [ {x[2]:.3f}]\n")