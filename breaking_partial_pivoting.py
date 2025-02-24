import all_gaussian_elimination_algorithms as GE
import numpy as np

# Here we devise a nonrandom matrix for which complete pivoting is significantly
# more accurate than partial pivoting.
n = 500
A_list = []
for i in range(n):
    row_i = [0] * n
    row_i[-1] = 1
    row_i[i] = 1
    for j in range(i):
        row_i[j] = -1
    A_list.append(row_i)

print(f"""
The matrix we are testing is of the form below, but of size {n}x{n}.

 [[  1,  0,  0,  0, 1],
  [ -1,  1,  0,  0, 1],
  [ -1, -1,  1,  0, 1],
  [ -1, -1, -1,  1, 1],
  [ -1, -1, -1, -1, 1]]
""")

print("We set b to be the product A @ x, where x = [1, 1, 1, ..., 1]")

A_numpy = np.array(A_list)

true_x = [1] * n
x_numpy = np.array(true_x)

b_numpy = A_numpy @ x_numpy
b_list = b_numpy.tolist()

partial_piv_x = GE.partial_piv_gaussian_elimination(A_list, b_list)
complete_piv_x = GE.complete_piv_gaussian_elimination(A_list, b_list)

err_partial_piv = float(np.linalg.norm(np.array(partial_piv_x) - x_numpy) / np.linalg.norm(x_numpy))
err_complete_piv = float(np.linalg.norm(np.array(complete_piv_x) - x_numpy) / np.linalg.norm(x_numpy))

print(f"Norm-wise relative error with partial pivoting: {err_partial_piv}.\n"
      f"Norm-wise relative error with complete pivoting: {err_complete_piv}.")