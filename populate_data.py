import random

# Specify the size of the square matrices
N = 1024 # You can change N to any size you want

# Generate random integers for matrices A and B
matrix_A = [[random.randint(0, 100) for _ in range(N)] for _ in range(N)]
matrix_B = [[random.randint(0, 100) for _ in range(N)] for _ in range(N)]

# Write the size N and the matrices to input00.txt
with open("input00.txt", "w") as file:
    # Write the size of the matrix N
    file.write(f"{N}\n")

    # Write matrix A
    for row in matrix_A:
        file.write(" ".join(map(str, row)) + "\n")

    # Write matrix B
    for row in matrix_B:
        file.write(" ".join(map(str, row)) + "\n")

print(f"File 'input00.txt' has been created with N = {N}")
