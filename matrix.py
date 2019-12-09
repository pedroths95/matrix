from random import seed
from random import randint
from random import random
import numpy
from numpy import copy
from numpy.linalg import inv
seed(23)
matrix=[]
n=7
for i in range(n):
    row=[]
    for j in range(n):
#        while rand_int != 0:
         f=random()
         row.append("{0:.2f}".format(f))
#        f=randint(0,n)
#        row.append('{0:+7.3f}'.format(f))
#        elif rand_int == 0:
    matrix.append(row)
#print(matrix)

def make_square_matrix(n): # n rows m colunms matrix
    matrix=[]
    for i in range(n):
        row=[]
        for j in range(n):
            f=randint(0,n)
            row.append(f)
        matrix.append(row)
    return matrix
def test_matrix(A):
    for i,l1 in enumerate(A):
        for j,l2 in enumerate(l1):
            if A[i][j] == A[j][i]:
                print ("Matrix element ",A[i][j],' is equal to',A[j][i],'for row = ',i,' and colunm = ',j)
            else:
                print("Matrix is Non Symmetric!")
# From https://github.com/ThomIves/BasicLinearAlgebraToolsPurePy/blob/master/LinearAlgebraPurePython.py
def zeros_matrix(rows, cols):
    """
    Creates a matrix filled with zeros.
        :param rows: the number of rows the matrix should have
        :param cols: the number of columns the matrix should have
        :return: list of lists that form the matrix
    """
    M = []
    while len(M) < rows:
        M.append([])
        while len(M[-1]) < cols:
            M[-1].append(0.0)

    return M


def copy_matrix(M):
    """
    Creates and returns a copy of a matrix.
        :param M: The matrix to be copied
        :return: A copy of the given matrix
    """
    # Section 1: Get matrix dimensions
    rows = len(M)
    cols = len(M[0])

    # Section 2: Create a new matrix of zeros
    MC = zeros_matrix(rows, cols)

    # Section 3: Copy values of M into the copy
    for i in range(rows):
        for j in range(cols):
            MC[i][j] = M[i][j]

    return MC
# https://integratedmlai.com/find-the-determinant-of-a-matrix-with-pure-python-without-numpy-or-scipy/
def determinant_recursive(A, total=0):
    # Section 1: store indices in list for row referencing
    indices = list(range(len(A)))
     
    # Section 2: when at 2x2 submatrices recursive calls end
    if len(A) == 2 and len(A[0]) == 2:
        val = A[0][0] * A[1][1] - A[1][0] * A[0][1]
        return val
 
    # Section 3: define submatrix for focus column and 
    #      call this function
    for fc in indices: # A) for each focus column, ...
        # find the submatrix ...
        As = copy_matrix(A) # B) make a copy, and ...
        As = As[1:] # ... C) remove the first row
        height = len(As) # D) 
 
        for i in range(height): 
            # E) for each remaining row of submatrix ...
            #     remove the focus column elements
            As[i] = As[i][0:fc] + As[i][fc+1:] 
 
        sign = (-1) ** (fc % 2) # F) 
        # G) pass submatrix recursively
        sub_det = determinant_recursive(As)
        # H) total all returns from recursion
        total += sign * A[0][fc] * sub_det 
 
    return total
# https://stackoverflow.com/questions/10508021/matrix-multiplication-in-pure-python
def matrixmult(A,B):
    rows_A = len(A)
    cols_A = len(A[0])
    rows_B = len(B)
    cols_B = len(B[0])

    if cols_A != rows_B:
      print ("Cannot multiply the two matrices. Incorrect dimensions.")
      return

    # Create the result matrix
    # Dimensions would be rows_A x cols_B
    C = [[0 for row in range(cols_B)] for col in range(rows_A)]
#    print (C)

    for i in range(rows_A):
        for j in range(cols_B):
            for k in range(cols_A):
                C[i][j] += A[i][k] * B[k][j]
    return C
def print_matrix(M,name):
    print(name)
    for row in M:
        print(row)
    print()
def check_symmetric(a, rtol=1e-05, atol=1e-08):
    return numpy.allclose(a, a.T, rtol=rtol, atol=atol)


A=make_square_matrix(4) # Our matrix is square
A=numpy.array(A)
print(check_symmetric(A))
if check_symmetric(A) == True:
    print('Symmetric Matrix!')
    exit()
else:
    pass
# https://www.geeksforgeeks.org/transpose-matrix-single-line-python/
A_t = [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]
print_matrix(A,'A')
print_matrix(A_t,'A_t')
letra_A=matrixmult(A,A_t)
Ai=inv(A)
letra_B=Ai
b=[randint(1,10) for i in range(4)]
print_matrix(b,'b')
letra_C=numpy.linalg.solve(A,b)
print_matrix(letra_A,'A*A^{t}')
print_matrix(letra_B,'A^{i}')
print_matrix(letra_C,'x')






print_matrix(matrixmult(A,Ai),'A*A_i')
print_matrix(matrixmult(Ai,A),'A_i*A')
EIG=numpy.linalg.eig
Ae=EIG(A)
print_matrix(Ae,' A Eigenvalues')
DET=numpy.linalg.det
print('Determinant calculated with numpy and pure python ',DET(A),determinant_recursive(A))
#print(seed())


