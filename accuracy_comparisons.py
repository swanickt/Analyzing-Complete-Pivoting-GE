import all_gaussian_elimination_algorithms as GE
import numpy as np

def random_accuracy_and_residual_checker(size: int) -> list[tuple[float, float]]:
    """
    Generate a random square coefficient matrix A of the given size. Set b such that
    the solution to Ax = b yields x_i = (-1)^{i+1}, where x = [x1, x2, ..., x_n]. Return
    a list of three tuples. The first holds the norm of the residual and norm-wise
    relative error for the solution obtained by Gaussian elimination with no pivoting.
    The second holds the same but for the solution obtained with partial pivoting, and the
    third for complete pivoting.
    """
    # Generate a random coefficient matrix (each coefficient is randomly chosen from [-30, 30])
    A = np.random.uniform(-30, 30, (size, size))


    # Compute b such that the solution x satisfies x_i = (-1)^{i+1} (where x = [x1, x2, ..., xn])
    # Note that we do (-1)^{i} below because of zero-based indexing
    true_x_value = np.array([(-1)**i for i in range(size)], dtype=float)
    b = A @ true_x_value

    # Cast the numpy arrays to nested lists
    A_list = A.tolist()
    b_list = b.tolist()

    # Use each algorithm to obtain a floating point approximation to true_x_value
    no_piv_x = GE.no_piv_gaussian_elimination(A_list, b_list)
    partial_piv_x = GE.partial_piv_gaussian_elimination(A_list, b_list)
    complete_piv_x = GE.complete_piv_gaussian_elimination(A_list, b_list)

    # Calculate the norm of the residuals for each solution
    r_no_piv = float(np.linalg.norm(b - (A @ np.array(no_piv_x))))
    r_partial_piv = float(np.linalg.norm(b - (A @ np.array(partial_piv_x))))
    r_complete_piv = float(np.linalg.norm(b - (A @ np.array(complete_piv_x))))

    # Calculate the accuracy of each solution, using their relative errors
    err_no_piv = float(np.linalg.norm(np.array(no_piv_x) - true_x_value) / np.linalg.norm(true_x_value))
    err_partial_piv = float(np.linalg.norm(np.array(partial_piv_x) - true_x_value) / np.linalg.norm(true_x_value))
    err_complete_piv = float(np.linalg.norm(np.array(complete_piv_x) - true_x_value) / np.linalg.norm(true_x_value))

    return [(r_no_piv, err_no_piv), (r_partial_piv, err_partial_piv), (r_complete_piv, err_complete_piv)]

print(" ----------------------------------------------------------------------------------------------------------------")
print("|System Size |No Piv Res. |Partial Piv Res. |Complete Piv Res. |No Piv Err. |Partial Piv Err. |Complete Piv Err. |")
print(" ----------------------------------------------------------------------------------------------------------------")

for n in range(5, 506, 20):
    results = random_accuracy_and_residual_checker(n)
    print("| {:<3}        | {:<3.4e} | {:<15.4e} | {:<16.4e} | {:<3.4e} | {:<15.4e} | {:<16.4e} |".
    format(n, results[0][0], results[1][0], results[2][0], results[0][1], results[1][1], results[2][1]))
print(" ----------------------------------------------------------------------------------------------------------------")